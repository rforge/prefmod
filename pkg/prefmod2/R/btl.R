## simple formula interface
btl <- function(formula, data, subset, na.action, weights, offset,
  type = c("loglin", "logit", "probit"), ref = NULL, ties = NULL,
  model = FALSE, y = TRUE, x = FALSE, ...)
{
  ## remember original call
  cl <- match.call()

  ## formula parsing
  m <- match.call(expand.dots = FALSE)
  m <- m[c(1, match(c("formula", "data", "subset", "weights", "na.action"), names(m), 0))]
  m[[1]] <- as.name("model.frame")
  mf <- eval.parent(m)

  ## extract data
  yy <- model.response(mf)  
  ww <- model.weights(mf)
  xx <- NULL ## FIXME: covariates not yet implemented

  ## call workhorse
  rval <- btl.fit(x = xx, y = yy, weights = ww,
    type = type, ref = ref, ties = ties, ...)

  ## add usual fitted model properties
  rval$call <- cl
  rval$terms <- terms(formula, data = data)
  rval$y <- if(y) yy else NULL
  rval$x <- if(x) xx else NULL
  rval$model <- if(model) mf else NULL

  return(rval)
}

## workhorse function
btl.fit <- function(x, y, weights = NULL,
  type = c("loglin", "logit", "probit"), ref = NULL, ties = NULL,
  ...)
{
  ## main arguments
  if(missing(x)) x <- NULL
  if(!is.null(x)) stop("regressors not yet implemented")
  if(missing(y)) stop("response missing")
  stopifnot(inherits(y, "paircomp"))  

  ## basic paircomp properties
  lab <- labels(y)
  nsubj <- length(y)
  nobj <- length(lab)
  npc <- nobj * (nobj - 1)/2
  mscale <- mscale(y)
  if(max(abs(mscale)) > 1) stop("comparisons on likert scales not yet implemented")
  has_ties <- mscale[2] == 0
  ix <- which(upper.tri(diag(nobj)), arr.ind = TRUE)
  
  ## weights
  if(isTRUE(all.equal(weights, rep(1, nsubj)))) weights <- NULL
  if(!is.null(weights)) {
    stopifnot(length(weights) == nsubj)
    if(!isTRUE(all.equal(weights, round(weights)))) stop("only case weights (integer) allowed")
    weights <- as.integer(round(weights))
    yorig <- y
    y <- rep(y, weights)
    nsubj <- length(y)
  }

  ## further arguments
  type <- match.arg(type, c("loglin", "logit", "probit"))
  if(is.null(ref)) ref <- nobj
  if(is.character(ref)) ref <- match(ref, lab)
  if(is.null(ties)) ties <- has_ties
  if(ties & type != "loglin") stop("only log-linear model can handle ties")  
  if(!has_ties & ties) stop("data have no ties")
  npar <- nobj - !ties
  
  ## basic aggregation quantities
  ytab <- summary(y)
  if(has_ties) ytab <- ytab[, c(1, 3, 2)]
  if(!ties) ytab <- ytab[, 1:2]
    
  ## set up auxiliary model
  if(type == "loglin") {
    famaux <- poisson()
    yaux <- as.vector(t(ytab))
    xaux <- matrix(0, nrow = 3 * npc, ncol = nobj + npc)
    for(i in 1:nrow(ix)) {
      ## xaux[i*3 - (2:1), ix[i,1:2]] <- c(1, -1, -1, 1) ## DHK parametrization
      xaux[i*3 - (2:0), ix[i,1:2]] <- c(1, 0, 0.5, 0, 1, 0.5)
      xaux[i*3 - (2:0), nobj + i] <- 1
    }
    xaux[,1:nobj] <- xaux[,c((1:nobj)[-ref], ref)]    
    xaux[,nobj] <- rep(c(0, 0, 1), npc)
    if(!ties) {
      xaux <- xaux[,-nobj]
      xaux <- xaux[-((1:npc) * 3),]
    }
  } else {
    famaux <- binomial(link = type)
    yaux <- ytab
    xaux <- matrix(0, nrow = npc, ncol = nobj)
    for(i in 1:nrow(ix)) xaux[i, ix[i,1:2]] <- c(1, -1)
    xaux <- xaux[,-ref]    
  }

  ## fit auxiliary model and extract information
  fm <- glm.fit(xaux, yaux, family = famaux, control = glm.control(...))
  par <- fm$coefficients[1:npar]
  vc <- summary.glm(fm, corr = FALSE)$cov.unscaled[1:npar,1:npar]
  names(par) <- rownames(vc) <- colnames(vc) <- c(lab[-ref], if(ties) "(undecided)" else NULL)

  ## log-probabilities and log-likelihood
  par2logprob <- switch(type,
    "loglin" = function(i) {
      p <- rep(0, nobj)
      p[-ref] <- if(ties) par[-npar] else par
      p <- p[ix[i,]]
      if(ties) p <- c(p, par[npar] + mean(p))
      p - log(sum(exp(p)))
    },
    "logit" = function(i) {
      p <- rep(0, nobj)
      p[-ref] <- par
      plogis(c(-1, 1) * diff(p[ix[i,]]), log.p = TRUE)
    },
    "probit" = function(i) {
      p <- rep(0, nobj)
      p[-ref] <- par
      pnorm(c(-1, 1) * diff(p[ix[i,]]), log.p = TRUE)
    }
  )
  logp <- t(sapply(1:npc, par2logprob))
  loglik <- sum(logp * ytab)

  ## raw original data
  if(!is.null(weights)) y <- yorig
  ymat <- as.matrix(y)

  ## estimating functions
  if(type == "probit") {
    ## not yet implemented
    ef <- NULL
  } else {
    if(!ties) logp <- cbind(logp, -Inf) ## ties impossible
    gradp <- matrix(0, nrow = npc * 3, ncol = nobj)
    cf <- -matrix(c(1, 0, 0, 0, 1, 0, 0.5, 0.5, 1), ncol = 3)
    ct <- matrix(c(1, 0, 0.5, 0, 1, 0.5, 0, 0, 1), ncol = 3)
    for(i in 1:npc) gradp[i*3 - (2:0), c(ix[i,], nobj)] <- t(t(ct) + as.vector(cf %*% exp(logp[i,])))
    ef <- t(sapply(1:length(y), function(i) {
      wi <- (0:(npc - 1)) * 3 + c(2, 3, 1)[ymat[i,] + 2]
      colSums(gradp[wi,], na.rm = TRUE)
    }))
    if(!ties) ef <- ef[,-nobj]
    dimnames(ef) <- list(names(y), names(par))
  }
  if(!is.null(weights)) ef <- ef * weights  
  
  ## collect, class, and return
  rval <- list(
    coefficients = par,
    vcov = vc,
    loglik = loglik,
    df = npar,
    estfun = ef,
    weights = weights,
    n = nsubj,
    type = type,
    ref = lab[ref],
    ties = ties,
    labels = lab
  )     
  class(rval) <- "btl"
  return(rval)
}

## standard methods
vcov.btl <- function(object, ...) object$vcov

logLik.btl <- function(object, ...) structure(object$loglik, df = object$df, class = "logLik")

deviance.btl <- function(object, ...) -2 * object$loglik

estfun.btl <- function(x, ...) {
  ef <- x$estfun
  if(is.null(ef)) stop(sprintf("estimating functions not yet implemented for %s BTL models"), x$type)
  return(ef)
}

bread.btl <- function(x, ...) x$vcov * x$n

## more elaborate methods
## FIXME: first draft only!
print.btl <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  if(is.null(x$call)) {
    cat("\nBradley-Terry-Luce model\n\n")  
  } else {
    cat("\nCall:\n")
    cat(paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  }
  cat("Parameters:\n")
  print(coef(x), digits = digits)
  cat("\nAssociated worth parameters:\n")
  print(worth(x), digits = digits)
  cat("\n")
  invisible(x)
}

summary.btl <- function(object, vcov. = NULL, ...)
{
  ## coefficients
  cf <- coef(object)

  ## covariance matrix
  if(is.null(vcov.)) 
      vc <- vcov(object)
  else {
      if(is.function(vcov.)) vc <- vcov.(object)
        else vc <- vcov.
  }
  
  ## Wald test of each coefficient
  cf <- cbind(cf, sqrt(diag(vc)), cf/sqrt(diag(vc)), 2 * pnorm(-abs(cf/sqrt(diag(vc)))))
  colnames(cf) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")

  object$coefficients <- cf      
  class(object) <- "summary.btl"
  return(object)
}

print.summary.btl <- function(x, digits = max(3, getOption("digits") - 3), 
    signif.stars = getOption("show.signif.stars"), ...)
{
  if(is.null(x$call)) {
    cat("\nBradley-Terry-Luce model\n\n")  
  } else {
    cat("\nCall:\n")
    cat(paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  }

  cat("Parameters:\n")
  printCoefmat(x$coefficients, digits = digits, signif.stars = signif.stars, na.print = "NA", ...)

  cat("\nLog-likelihood:", format(signif(x$loglik, digits)),
    "(df =", paste(x$df, ")", sep = ""), "\n\n")
  invisible(x)
}


coef.btl <- function(object, all = TRUE, ref = !all, ...) {
  lab <- object$labels
  nobj <- length(lab)
  acf <- object$coefficients
  cf <- structure(rep(0, nobj), .Names = lab)
  ocf <- acf[1:(nobj-1)]
  cf[names(ocf)] <- ocf
  cf <- c(cf, acf[-(1:(nobj-1))])
  if(!all) cf <- cf[1:nobj]
  if(!ref) cf <- cf[-match(object$ref, lab)]
  return(cf)
}

worth <- function(object, ...) UseMethod("worth")

worth.btl <- function(x, ...) {
  if(x$type == "probit") stop("Worth parameters for probit model not (yet) available.")
  lab <- x$labels
  cf <- coef(x, all = FALSE, ref = TRUE)
  exp(cf)/sum(exp(cf))
}

plot.btl <- function(x, 
  worth = TRUE, index = TRUE, names = TRUE, ref = TRUE, abbreviate = FALSE,
  type = NULL, lty = NULL, xlab = "Objects", ylab = NULL, ...)
{
  ## parameters to be plotted
  if(worth) {
    cf <- worth(x)
  } else {
    cf <- coef(x, all = FALSE, ref = TRUE)
    if(is.character(ref) | is.numeric(ref)) {
      reflab <- ref
      ref <- TRUE
    } else {
      reflab <- x$ref
    }
    if(is.character(reflab)) reflab <- match(reflab, x$labels)
    cf <- cf - cf[reflab]
  }

  ## labeling
  cf_ref <- if(!worth) 0 else 1/length(cf)
  if(is.character(names)) {
    names(cf) <- names
    names <- TRUE
  }
  if(is.null(ylab)) ylab <- if(worth) "Worth parameters" else "Parameters"
    
  ## abbreviation
  if(is.logical(abbreviate)) {
    nlab <- max(nchar(names(cf)))
    abbreviate <- if(abbreviate) as.numeric(cut(nlab, c(-Inf, 1.5, 4.5, 7.5, Inf))) else nlab
  }
  names(cf) <- abbreviate(names(cf), abbreviate)

  ## raw plot
  ix <- if(index) seq(along = cf) else rep(0, length(cf))
  plot(ix, cf, xlab = xlab, ylab = ylab, type = "n", axes = FALSE, ...)
  if(ref) abline(h = cf_ref, col = "lightgray")
  axis(2)
  box()  

  ## actual data
  if(index) {
    if(is.null(type)) type <- "b"
    if(is.null(lty)) lty <- 2  
    lines(ix, cf, type = type, lty = lty)
    axis(1, at = ix, labels = if(names) names(cf) else TRUE)
  } else {
    if(is.null(type)) type <- "p"
    if(names) text(names(cf), x = ix, y = cf, ...) else lines(ix, cf, type = type, lty = lty)
  }
}
