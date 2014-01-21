#############################################################################################################################
## Temporary copies of several internal psychotools functions.
## To be deleted with the next psychotools release
#############################################################################################################################


itempar.RaschModel <- function (object, ref = NULL, vcov = TRUE, ...) {

  ## extract cf and labels, include restricted parameter
  cf <- c(0.00, coef(object))
  m <- length(cf)
  lbs <- names(cf)
  lbs[1] <- ifelse(lbs[1] == "Item02", "Item01", colnames(object$data)[1])

  ## process ref
  if (is.null(ref)) {
    ref <- 1:m
  } else if (is.character(ref)) {
    stopifnot(all(ref %in% lbs))
    ref <- which(lbs %in% ref)
  } else if (is.numeric(ref)) {
    ref <- as.integer(ref)
    stopifnot(all(ref %in% 1:m))
  } else stop("Argument 'ref' can only be a character vector with item labels or a numeric vector with item indices.")

  ## create contrast matrix
  D <- diag(m)
  D[, ref] <- D[, ref] - 1/length(ref)

  ## impose restriction
  cf <- as.vector(D %*% cf)
  names(cf) <- lbs

  ## create adjusted vcov if requested
  if (vcov) {
    vc <- D %*% rbind(0, cbind(0, vcov(object))) %*% t(D)
    colnames(vc) <- rownames(vc) <- lbs
  }

  ## return results
  rval <- structure(cf, class = "itempar", model = "RM", ref = ref, vcov = if (vcov) vc else NULL)
  return(rval)
}

itempar.RSModel <- function (object, ref = NULL, vcov = TRUE, simplify = TRUE, ...) {

  ## extract cf and labels, include restricted parameters
  cf <- coef(object)
  k <- length(cf)
  m <- sum(object$items)
  icf <- c(0.00, cf[1:(m-1)])
  ccf <- c("C1" = 0.00, cf[m:k])
  ilbs <- names(icf)
  clbs <- names(ccf)
  names(cf)[1] <- names(icf)[1] <- ilbs[1] <- if (ilbs[2] == "I2") "I1" else colnames(object$data)[1]
  
  ## process ref
  if (is.null(ref)) {
    ref <- 1:m
  } else if (is.character(ref)) {
    stopifnot(all(ref %in% ilbs))
    ref <- which(ilbs %in% ref)
  } else if (is.numeric(ref)) {
    ref <- as.integer(ref)
    stopifnot(all(ref %in% 1:m))
  } else stop("Argument 'ref' can only be a character vector with item labels or a numeric vector with item indices.")

  ## create contrast matrix
  D <- diag(m)
  D[, ref] <- D[, ref] - 1/length(ref)
  
  ## impose restriction
  icf <- D %*% icf
  names(icf) <- ilbs

  ## adjust vcov if requested
  if (vcov) {
    k <- k + 2
    vc <- matrix(0.00, nrow = k, ncol = k)
    vc[c(2:m, (m+2):k), c(2:m, (m+2):k)] <- vcov(object)

    ## create new D to include covariances of ip and taus in transformation
    D2 <- diag(k)                       
    D2[1:m, 1:m] <- D
    vc <- D2 %*% vc %*% t(D2)
    colnames(vc) <- rownames(vc) <- c(ilbs, clbs)
  }

  ## create list if requested
  cf <- if (simplify) c(icf, ccf) else list(icf, ccf)

  ## return results
  rval <- structure(cf, class = "itempar", model = "RSM", ref = ref, vcov = if (vcov) vc else NULL)
  return(rval)
}

itempar.PCModel <- function (object, ref = NULL, vcov = TRUE, simplify = TRUE, ...) {
  
  ## extract cf and labels, include restricted parameter
  cf <- c(0.00, coef(object))
  m <- length(cf)
  lbs <- names(cf)
  lbs[1] <- names(cf)[1] <- if (lbs[2] == "I1-C2") "I1-C1" else if (lbs[2] == "I2-C1") "I1-C1" else paste0(colnames(object$data)[1], "-C1")
  
  ## process ref
  if (is.null(ref)) {
    ref <- 1:m
  } else if (is.character(ref)) {
    stopifnot(all(ref %in% lbs))
    ref <- which(lbs %in% ref)
  } else if (is.numeric(ref)) {
    ref <- as.integer(ref)
    stopifnot(all(ref %in% 1:m))
  } else stop("Argument 'ref' can only be a character vector with threshold parameter labels or a numeric vector with threshold parameter indices.")

  ## create contrast matrix
  D <- diag(m)
  D[, ref] <- D[, ref] - 1/length(ref)

  ## impose restriction
  cf <- D %*% cf
  names(cf) <- lbs

  ## create adjusted vcov if requested
  if (vcov) {
    vc <- D %*% rbind(0, cbind(0, vcov(object))) %*% t(D)
    colnames(vc) <- rownames(vc) <- lbs
  }

  ## create list if requested
  if (!simplify) {
    cf <- relist(cf, object$categories)
    lbs <- relist(lbs, object$categories)
    cf <- mapply("names<-", cf, lbs, SIMPLIFY = FALSE, USE.NAMES = FALSE)
  }

  ## return results
  rval <- structure(cf, class = "itempar", model = "PCM", ref = ref, vcov = if (vcov) vc else NULL)
  return(rval)
}

threshold.RaschModel <- function (object, type = c("mode", "median", "mean", "unmodified"), ref = NULL, ...) {
  ## check type, extract item parameters
  type <- match.arg(type)
  ip <- itempar.RaschModel(object, ref = ref, vcov = FALSE)

  ## return requested threshold parameters // in RM: mode = median = mean = unmodified
  class(ip) <- "threshold"
  attr(ip, "type") <- type
  return(ip)
}

threshold.RSModel <- function (object, type = c("mode", "median", "mean", "unmodified"), ref = NULL, simplify = TRUE, ...) {

  ## check type, extract item parameters, add names, backup attributes
  type <- match.arg(type)
  ip <- itempar.RSModel(object, ref = ref, vcov = FALSE, simplify = FALSE)
  lbs <- lapply(names(ip[[1]]), paste, names(ip[[2]]), sep = "-")
  mdl <- attr(ip, "model")
  ref <- attr(ip, "ref")
  ip <- lapply(as.list(ip[[1]]), function (beta) diff(0:length(ip[[2]]) * beta + c(0, ip[[2]]))) ## convert to pcm parameters, rest as before
  ip <- mapply("names<-", ip, lbs, SIMPLIFY = FALSE)
  names(ip) <- NULL

  ## calculate threshold parameters (if necessary)
  if (type == "mode") {
   ## check if all threshold parameters are in order, if not, calculate sorted ones
    us <- sapply(ip, is.unsorted)
    if (any(us)) {
      usj <- which(us)
      for (j in usj) {
        ipj <- ip[[j]]
        nj <- length(ipj)
        
        ## check if there is a point with a biggest parameter, if yes, take mean
        for (i in 1:nj) {
          if (all(ipj[i] > ipj[(i+1):nj])) {
            ipj[i] <- mean(ipj[i:nj])
            ipj <- ipj[-(i+1:nj)]
            break
          }
        }
        
        ## recursive sorting if there is still unorder (e.g. 4, 2, 3, 1)
        while(is.unsorted(ipj)) {
          uo_pos <- which(diff(ipj) < 0)                             # locate unordered parameters, returns position of the first
          ipj[uo_pos] <- (ipj[uo_pos] + ipj[uo_pos + 1]) / 2 # replace first with location of intersection of ccc curves (= (eps1 + eps2)/ 2)
          ipj <- ipj[-(uo_pos + 1)]                              # remove second
        }

        ip[[j]] <- ipj
      }
    }
  } else {
    oj <- sapply(ip, length)

    if (type == "median") {
      ## function to find locations on theta axis
      zmedian <- function (theta = NULL, delta = NULL, geq = NULL, ncat = NULL) {
        rowSums(ppcm(theta = theta, delta = delta)[, (geq + 1):ncat, drop = FALSE]) - 0.5
      }
      
      ## loop though items and find locations by means of zmedian() and uniroot()
      for (j in seq_along(ip)) {
        ip[[j]] <- sapply(1:oj[j], function (geq) uniroot(f = zmedian, interval = c(-10, 10), delta = ip[[j]], geq = geq, ncat = oj[j] + 1)$root)
      }
    } else if (type == "mean") {
      ## function to find locations on theta axis
      xpct <- lapply(oj, function (oj) 1:oj - 0.5)
      zexpct <- function (theta = NULL, delta = NULL, expct = NULL) ppcm(theta = theta, delta = delta) %*% 0:length(delta) - expct
      
      ## loop though items and find locations by means of zexpct() and uniroot()
      for (j in seq_along(ip)) {
        ip[[j]] <- sapply(xpct[[j]], function (xp) uniroot(f = zexpct, interval = c(-10, 10), delta = ip[[j]], expct = xp)$root)
      }
    }

    ## add labels again
    ip <- mapply("names<-", ip, lbs, SIMPLIFY = FALSE)
  }
  
  ## backup attributes, simplify if requested, then set attributes again
  if (simplify) ip <- unlist(ip)
  return(structure(ip, "class" = "threshold", "model" = mdl, "ref" = ref, "type" = type))
}

threshold.PCModel <- function (object, type = c("mode", "median", "mean", "unmodified"), ref = NULL, simplify = TRUE, ...) {
  ## check type, extract item parameters, backup attributes
  type <- match.arg(type)
  ip <- itempar.PCModel(object, ref = ref, vcov = FALSE, simplify = FALSE)
  mdl <- attr(ip, "model")
  ref <- attr(ip, "ref")
  ip <- lapply(ip, function (j) diff.default(c(0, j))) # item category parameters -> item threshold parameters

  ## calculate threshold parameters (if necessary)
  if (type == "mode") {
    ## check if all threshold parameters are in order, if not, calculate sorted ones
    us <- sapply(ip, is.unsorted)
    if (any(us)) {
      usj <- which(us)
      for (j in usj) {
        ipj <- ip[[j]]
        nj <- length(ipj)
        
        ## check if there is a point with a biggest parameter, if yes, take mean
        for (i in 1:nj) {
          if (all(ipj[i] > ipj[(i+1):nj])) {
            ipj[i] <- mean(ipj[i:nj])
            ipj <- ipj[-(i+1:nj)]
            break
          }
        }
        
        ## recursive sorting if there is still unorder (e.g. 4, 2, 3, 1)
        while(is.unsorted(ipj)) {
          uo_pos <- which(diff(ipj) < 0)                             # locate unordered parameters, returns position of the first
          ipj[uo_pos] <- (ipj[uo_pos] + ipj[uo_pos + 1]) / 2 # replace first with location of intersection of ccc curves (= (eps1 + eps2)/ 2)
          ipj <- ipj[-(uo_pos + 1)]                              # remove second
        }

        ip[[j]] <- ipj
      }
    }
  } else {
    ## backup labels
    lbs <- lapply(ip, names)
    oj <- sapply(ip, length)

    if (type == "median") {
      ## function to find locations on theta axis
      zmedian <- function (theta = NULL, delta = NULL, geq = NULL, ncat = NULL) {
        rowSums(ppcm(theta = theta, delta = delta)[, (geq + 1):ncat, drop = FALSE]) - 0.5
      }
      
      ## loop though items and find locations by means of zmedian() and uniroot()
      for (j in seq_along(ip)) {
        ip[[j]] <- sapply(1:oj[j], function (geq) uniroot(f = zmedian, interval = c(-10, 10), delta = ip[[j]], geq = geq, ncat = oj[j] + 1)$root)
      }
    } else if (type == "mean") {
      ## function to find locations on theta axis
      xpct <- lapply(oj, function (oj) 1:oj - 0.5)
      zexpct <- function (theta = NULL, delta = NULL, expct = NULL) ppcm(theta = theta, delta = delta) %*% 0:length(delta) - expct
      
      ## loop though items and find locations by means of zexpct() and uniroot()
      for (j in seq_along(ip)) {
        ip[[j]] <- sapply(xpct[[j]], function (xp) uniroot(f = zexpct, interval = c(-10, 10), delta = ip[[j]], expct = xp)$root)
      }
    }

    ## add labels again
    ip <- mapply("names<-", ip, lbs, SIMPLIFY = FALSE)
  }

  ## backup attributes, simplify if requested, then set attributes again
  if (simplify) ip <- unlist(ip)
  return(structure(ip, "class" = "threshold", "model" = mdl, "ref" = ref, "type" = type))
}

ppcm <- function(theta = NULL, delta = NULL)
{
  ## check input
  stopifnot(!is.null(theta) && !is.null(delta))
  
  ## if list input, recurse...
  if (is.list(theta)) return(lapply(theta, ppcm, delta = delta))
  if (is.list(delta)) return(lapply(delta, ppcm, theta = theta))

  ## calculate probabilities
  num <- cbind(0, outer(theta, delta, "-")) # all possible differences, 0 for category zero (\sum_0^0 \def 0)
  num <- t(exp(apply(num, 1, cumsum)))      # numerator: all possible cumulative sums
  denom <- rowSums(num)                     # denominator: sum over cumulative sums
  return(num/denom)
}
