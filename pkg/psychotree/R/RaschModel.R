## S4 StatModel object
RaschModel <- function(gradtol = 1e-12, hessian = FALSE) {
  new("StatModel",
    capabilities = new("StatModelCapabilities"),
    name = "Rasch model",
    dpp = ModelEnvFormula,
    fit = function(object, weights = NULL, ...){

        ## extract response (there are no regressors)
        y <- object@get("response")

        ## call RaschModel.fit()
        z <- RaschModel.fit(y = y, weights = weights,
          gradtol = gradtol, hessian = hessian)
        z$ModelEnv <- object
        z$addargs <- list(...)
        z
    }
  )
}

## methods needed for mob()
reweight.RaschModel <- function(object, weights, ...) {
     fit <- RaschModel(gradtol = object$gradtol, hessian = object$hessian)@fit
     do.call("fit", c(list(object = object$ModelEnv, weights = weights), object$addargs))
}

print.RaschModel <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("Rasch model coefficients:\n")
  print(coef(x), digits = digits)
  invisible(x)
}

## workhorse fitting function
RaschModel.fit <- function(y, weights = NULL, start = NULL,
  gradtol = 1e-12, check.analyticals = FALSE, hessian = FALSE, ...)
{
  ## original data
  y <- as.matrix(y)
  k <- k_orig <- ncol(y)
  n <- nrow(y)
  if(is.null(colnames(y))) colnames(y) <- paste("Item", gsub(" ", "0", format(1:k)), sep = "")
  y_orig <- y

  ## weights processing
  if(is.null(weights)) weights <- rep.int(1L, n)
  ## data and weights need to match
  stopifnot(length(weights) == n)

  ## all parameters identified?
  cm <- colMeans(y, na.rm = TRUE)
  status <- as.character(cut(cm, c(-Inf, 1/(2 * n), 1 - 1/(2 * n), Inf), labels = c("0", "0/1", "1")))
  status[is.na(status)] <- "NA"
  status <- factor(status, levels = c("0/1", "0", "1", "NA"))
  ident <- status == "0/1"
  names(status) <- colnames(y)

  ## just estimate identified parameters
  y <- y[,ident, drop = FALSE]
  k <- ncol(y)
  y_na <- is.na(y)
  any_y_na <- any(y_na)

  if(!any_y_na) {
    ## compute likelihood/gradient/hessian on aggregated data
  
    ## data statistics
    cs <- colSums(y * weights)
    rs <- rowSums(y)
    rf <- as.vector(tapply(weights, factor(rs, levels = 0:k), sum))
    rf[is.na(rf)] <- 0

    ## starting values
    if(is.null(start)) start <- qlogis(cs/sum(weights))
  
    ## contrast: set parameter 1 to zero
    start <- start[-1] - start[1]
    rf <- rf[-1]
    cs <- cs[-1]

    ## objective function: conditional log-likelihood
    cloglik <- function(par) {
      ## obtain esf and apply contrast
      esf <- elementary_symmetric_functions(c(0, par), order = 1)
      g <- esf[[1]][-1]
      g1 <- esf[[2]][-1, -1]

      ## conditional log-likelihood
      cll <- sum(cs * par) - sum(rf * log(g))
      ## gradient
      grad <- cs - colSums(rf * g1/g)

      ## collect and return
      attr(cll, "gradient") <- - grad
      return(-cll)
    }

    ## analytical hessian
    ahessian <- function(par, esf) {
      ## obtain esf and apply contrast
      g <- esf[[1]][-1]
      g1 <- esf[[2]][-1, -1]
      g2 <- esf[[3]][-1, -1, -1]

      ## hessian
      hess <- matrix(0, ncol = k-1, nrow = k-1)
      g1s <- g1/g
      for (q in 1:(k-1)) hess[q,] <- colSums(rf * (g2[,q,]/g - (g1[,q]/g) * g1s))
    
      return(hess)
    }

  } else {
    ## compute likelihood/gradient/hessian on individual data

    na_patterns <- factor(apply(y_na, 1, function(z) paste(which(z), collapse = "\r")))
    y0 <- y
    y0[y_na] <- 0
    rs <- rowSums(y0)

    ## starting values
    if(is.null(start)) start <- qlogis(colSums(y0 * weights)/colSums(!y_na * weights))
    start <- start[-1] - start[1]

    ## objective function: conditional log-likelihood
    cloglik <- function(par) {
      ## pre-compute ESF for all observed NA patterns
      esf_patterns <- lapply(levels(na_patterns), function(z) {
        wi <- as.integer(strsplit(z, "\r")[[1]])
	cf <- if(length(wi) < 1) c(0, par) else c(0, par)[-wi]
	elementary_symmetric_functions(cf, order = 1)
      })
      
      ## convenience function for obtaining ESF of observation i
      ## fill up with zeros if there were NAs
      get_esf <- function(i, order = 0) {
        na_i <- na_patterns[i]
        rs_i <- rs[i]
        wi_i <- as.integer(strsplit(as.character(na_i), "\r")[[1]])
        esf <- esf_patterns[[na_i]][[order + 1]]
	if(order < 1) esf[-1][rs_i] else {
	  rval <- rep(0, k)
	  if(length(wi_i) < 1) {
	    rval <- esf[rs_i + 1,]
	  } else {
            rval[-wi_i] <- esf[rs_i + 1,]
	  }
          rval[-1]
	}
      }

      ## compute sum of weighted likelihood contributions
      cll <- sum(weights * (as.vector(y0[,-1] %*% par) -
        sapply(1:n, function(i) log(get_esf(i, order = 0)))))

      ## gradient
      grad <- colSums(weights * (y0[,-1] - t(sapply(1:n, function(i)
        get_esf(i, order = 1) / get_esf(i, order = 0)))))

      ## collect and return
      attr(cll, "gradient") <- - grad
      return(-cll)
    }
 
    ## analytical gradient
    agrad <- function(par, esf_patterns) {
      ## convenience function for obtaining ESF of observation i
      ## fill up with zeros if there were NAs
      get_esf <- function(i, order = 0) {
        na_i <- na_patterns[i]
        rs_i <- rs[i]
        wi_i <- as.integer(strsplit(as.character(na_i), "\r")[[1]])
        esf <- esf_patterns[[na_i]][[order + 1]]
	if(order < 1) esf[-1][rs_i] else {
	  rval <- rep(0, k)
	  if(length(wi_i) < 1) {
	    rval <- esf[rs_i + 1,]
	  } else {
            rval[-wi_i] <- esf[rs_i + 1,]
	  }
          rval[-1]
	}
      }

      ## gradient
      weights * (y0[,-1] - t(sapply(1:n, function(i)
        get_esf(i, order = 1) / get_esf(i, order = 0))))
    }

    ## analytical hessian
    ahessian <- function(par, esf) {
      ## convenience function for obtaining ESF of observation i
      ## fill up with zeros if there were NAs
      esf_patterns <- esf
      get_esf <- function(i, order = 0) {
        na_i <- na_patterns[i]
        rs_i <- rs[i]
        wi_i <- as.integer(strsplit(as.character(na_i), "\r")[[1]])
        esf <- esf_patterns[[na_i]][[order + 1]]
	if(order < 1) {
	  esf[-1][rs_i]
	} else if(order < 2) {
	  rval <- rep(0, k)
	  if(length(wi_i) < 1) {
	    rval <- esf[rs_i + 1,]
	  } else {
            rval[-wi_i] <- esf[rs_i + 1,]
	  }
          rval[-1]
	} else {
	  rval <- matrix(0, ncol = k, nrow = k)
	  if(length(wi_i) < 1) {
	    rval <- esf[rs_i + 1,,]
	  } else {
            rval[-wi_i, -wi_i] <- esf[rs_i + 1,,]
	  }
          rval[-1,-1]
	}
      }

      hess <- matrix(0, ncol = k-1, nrow = k-1)

      for(i in 1:n) {
        g <- get_esf(i, order = 0)
        g1 <- get_esf(i, order = 1)
        g2 <- get_esf(i, order = 2)
        hess <- hess + (g2 - outer(g1, g1)/g)/g
      }

      return(hess)
    }

  }
  
  ## optimization
  opt <- nlm(cloglik, start,
    hessian = hessian, gradtol = gradtol, check.analyticals = check.analyticals)
  
  ## collect and annotate results
  cf <- opt$estimate
  esf <- if(any_y_na) {
    lapply(levels(na_patterns), function(z) {
      wi <- as.integer(strsplit(z, "\r")[[1]])
      cfi <- if(length(wi) < 1) c(0, cf) else c(0, cf)[-wi]
      elementary_symmetric_functions(cfi, order = 2)
    })
  } else {
    elementary_symmetric_functions(c(0, cf), order = 2)
  }
  grad <- if(any_y_na) agrad(cf, esf) else NULL
  vc <- if(hessian) opt$hessian else ahessian(cf, esf)
  vc <- solve(vc)
  names(cf) <- rownames(vc) <- colnames(vc) <- colnames(y)[-1]
  
  ## collect, class, and return
  rval <- list(
    coefficients = cf,
    vcov = vc,
    loglik = -opt$minimum,
    df = k-1,
    weights = weights,
    n = n,
    data = y_orig,
    items = status,
    na = any_y_na,
    elementary_symmetric_functions = if(any_y_na) NULL else esf,
    estfun = grad,
    nlm_code = opt$code,
    iterations = opt$iterations,
    hessian = hessian,
    gradtol = gradtol
  )
  class(rval) <- "RaschModel"
  return(rval)
}

## elementary symmetric functions
elementary_symmetric_functions <- function(par, order = 0, log = TRUE) {
  ## Michelle Liou (1994). More on the Computation of Higher-Order
  ## Derivatives of the Elementary Symmetric Functions in the
  ## Rasch Model. Applied Psychological Measurement, 18(1), 53-62.

  ## derivatives up to order
  order <- round(order)[1]
  stopifnot(order %in% 0:2)
  rval <- list()[1:(order+1)]
  names(rval) <- 0:order

  ## transformations
  par <- as.vector(par)
  beta <- if(log) par else log(par)
  eps <- exp(beta)
  n <- length(eps)
  stopifnot(n > 2)
  
  ## Order: 0  
  ## initialization: gamma_1(eps_1) = eps_1
  gamma <- c(eps[1], numeric(n-1))

  ## recursion: Equation 3
  for(i in 2:n) gamma[1:i] <- c(eps[i] + gamma[1],
    eps[i] * gamma[(2:i) - 1] + gamma[2:i])
  
  ## gamma_0 = 1
  gamma <- c(1, gamma)

  ## return value
  rval[[1]] <- gamma
  if(order < 1) return(rval)

  ## Order: 1
  ## initialization: gamma_1^(j) = 1
  gamma1 <- matrix(0, nrow = n+1, ncol = n)
  gamma1[2,] <- 1
  ## recursion: Equation 4
  for(q in 3:(n+1)) gamma1[q,] <- gamma[q-1] - eps * gamma1[q-1,]
  ## if input on log scale: include inner derivative
  if(log) gamma1 <- t(t(gamma1) * eps)
  ## return value
  rval[[2]] <- gamma1
  if(order < 2) return(rval)

  ## Order: 2
  ## initialization: gamma_2^(i,j) = 1
  gamma2 <- array(0, c(n+1, n, n))
  gamma2[3,,] <- 1
  ## auxiliary variables
  eps_plus_eps <- outer(eps, eps, "+")
  eps_times_eps <- outer(eps, eps)
  ## recursion: Jansen's Equation (Table 1, Forward, Second-Order)
  for(q in 4:(n+1)) gamma2[q,,] <- gamma[q-2] -
    (eps_plus_eps * gamma2[q-1,,] + eps_times_eps * gamma2[q-2,,])
  ## if input on log scale: include inner derivative
  if(log) for(q in 1:(n+1)) gamma2[q,,] <- eps_times_eps * gamma2[q,,]
  ## each diagonal is simply first derivative
  for(q in 2:(n+1)) diag(gamma2[q,,]) <- gamma1[q,]
  ## return value
  rval[[3]] <- gamma2
  return(rval)
}

## methods
coef.RaschModel <- function(object, ...) object$coefficients

vcov.RaschModel <- function(object, ...) object$vcov

logLik.RaschModel <- function(object, ...) structure(object$loglik, df = object$df, class = "logLik")

bread.RaschModel <- function(x, ...) x$vcov * x$n

estfun.RaschModel <- function(x, ...) {
  if(!x$na) {
    g <- x$elementary_symmetric_functions[[1]][-1]
    g1 <- x$elementary_symmetric_functions[[2]][-1,-1]
    dat <- x$data[, x$items == "0/1", drop = FALSE]
    rs <- rowSums(dat)
    return(dat[,-1] - g1[rs,] / g[rs])
  } else {
    return(x$estfun)
  }
}

worth.RaschModel <- function(object, ...) {
  cf <- c(0, object$coefficients)
  cf <- cf - mean(cf)
  rval <- structure(rep(NA, ncol(object$data)), .Names = colnames(object$data))
  rval[object$items == "0/1"] <- cf
  rval[object$items == "0"] <- -Inf
  rval[object$items == "1"] <- Inf
  return(rval)
}

summary.RaschModel <- function(object, vcov. = NULL, ...)
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
  class(object) <- "summary.RaschModel"
  return(object)
}

print.summary.RaschModel <- function(x, digits = max(3, getOption("digits") - 3), 
    signif.stars = getOption("show.signif.stars"), ...)
{
  if(is.null(x$call)) {
    cat("\nRasch model\n\n")  
  } else {
    cat("\nCall:\n")
    cat(paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  }

  if(any(x$items != "0/1")) cat("Excluded items:",
    paste(colnames(x$data)[x$items != "0/1"], collapse = ", "), "\n\n")

  cat("Parameters:\n")
  printCoefmat(x$coefficients, digits = digits, signif.stars = signif.stars, na.print = "NA", ...)

  cat("\nLog-likelihood:", format(signif(x$loglik, digits)),
    "(df =", paste(x$df, ")", sep = ""), "\n")
  cat("Number of iterations in nlm optimization:", x$iterations, "\n\n")
  invisible(x)
}

plot.RaschModel <- function(x,
  center = TRUE, index = TRUE, names = TRUE, abbreviate = FALSE, ref = TRUE,
  bg = hcl(c(0, 0, 180, 270), c(0, 80, 80, 80), c(95, 50, 50, 50)), cex = 1,
  type = NULL, lty = NULL, ylim = NULL, xlab = "Items", ylab = NULL, ...)
{
  ## parameters to be plotted
  cf <- worth(x)
  cf_ident <- is.finite(cf) & !is.na(cf)
  if(!center) cf <- cf - (cf[cf_ident])[1]
  cf_ref <- mean(cf[cf_ident])

  ## background color
  col <- if(index) bg[1] else hcl(0, 0, 0)
  col <- rep(col, length(cf))
  col[cf <= -Inf] <- bg[2]
  col[cf >= Inf] <- bg[3]
  col[is.na(cf)] <- bg[4]

  ## substitute non-identified parameters with plottable values
  if(is.null(ylim)) ylim <- range(cf[cf_ident])
  ylim <- rep(ylim, length.out = 2)
  cf[is.na(cf)] <- cf_ref
  cf[cf <= -Inf] <- ylim[1]
  cf[cf >= Inf] <- ylim[2]

  ## labeling
  if(is.character(names)) {
    names(cf) <- names
    names <- TRUE
  }

  if(is.null(ylab)) ylab <- if(center) "Centered item parameters" else "Item parameters"

  ## abbreviation
  if(is.logical(abbreviate)) {
    nlab <- max(nchar(names(cf)))
    abbreviate <- if(abbreviate) as.numeric(cut(nlab, c(-Inf, 1.5, 4.5, 7.5, Inf))) else nlab
  }
  names(cf) <- abbreviate(names(cf), abbreviate)

  ## raw plot
  ix <- if(index) seq(along = cf) else rep(0, length(cf))
  plot(ix, cf, xlab = xlab, ylab = ylab, type = "n", axes = FALSE, ylim = ylim, ...)
  if(ref) abline(h = cf_ref, col = "lightgray")
  axis(2)
  box()  

  ## actual data
  if(index) {
    if(is.null(type)) type <- "b"
    if(is.null(lty)) lty <- 2
    if(type %in% c("b", "p")) points(ix, cf, pch = 19, col = col, cex = cex)
    lines(ix, cf, type = type, lty = lty, col = ifelse(cf_ident, "black", "transparent"), cex = cex)
    axis(1, at = ix, labels = if(names) names(cf) else TRUE)
  } else {
    if(is.null(type)) type <- "p"
    if(names) text(names(cf), x = ix, y = cf, col = col, ...) else lines(ix, cf, type = type, lty = lty, col = col)
  }
}

