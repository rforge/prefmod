## simple formula interface
btl <- function(formula, data, subset, na.action, weights, offset,
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
  rval <- btl.fit(x = xx, y = yy, weights = ww, ...)

  ## add usual fitted model properties
  rval$call <- cl
  rval$terms <- terms(formula, data = data)
  rval$y <- if(y) yy else NULL
  rval$x <- if(x) xx else NULL
  rval$model <- if(model) mf else NULL
  class(rval) <- "btl"

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
      plogis(c(1, -1) * diff(p[ix[i,]]))
    },
    "probit" = function(i) {
      p <- rep(0, nobj)
      p[-ref] <- par
      pnorm(c(1, -1) * diff(p[ix[i,]]))
    }
  )
  logp <- t(sapply(1:npc, par2logprob))
  loglik <- sum(logp * ytab)

  ## raw original data
  ymat <- if(is.null(weights)) as.matrix(y) else as.matrix(yorig)

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
  return(rval)
}

## FIXME: first draft only!
print.btl <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("Bradley-Terry-Luce model\n\nCoefficients:\n")
  print(coef(x), digits = digits)
  cat("\nAssociated worth parameters:\n")
  print(worth(x), digits = digits)
  invisible(x)
}

vcov.btl <- function(object, ...) object$vcov

logLik.btl <- function(object, ...) structure(object$loglik, df = object$df, class = "logLik")

deviance.btl <- function(object, ...) -2 * object$loglik

estfun.btl <- function(x, ...) {
  ef <- x$estfun
  if(is.null(ef)) stop(sprintf("estimating functions not yet implemented for %s BTL models"), x$type)
  return(ef)
}

bread.btl <- function(x, ...) x$vcov * x$n

worth <- function(object, ...) UseMethod("worth")

worth.btl <- function(x, ...) {
  lab <- x$labels
  rval <- structure(rep(0, length(lab)), .Names = lab)
  cf <- coef(x)[1:(length(lab)-1)]
  rval[names(cf)] <- cf
  exp(rval)/sum(exp(rval))
}
