## paircomp() is the class constructor:
##   - data: integer coded matrix of n subjects (rows) and
##            choose(k, 2) paired comparisons (columns)
##   - labels: labels of k objects
##   - mscale: measurement scale for comparisons (centered around 0)
##   - ordered: logical. Is a:b different from b:a?
##   - covariates: data frame with k rows for object covariates

paircomp <- function(data, labels = NULL, mscale = NULL, ordered = FALSE, covariates = NULL)
{
  ## data should be matrix(-like)
  npc <- NCOL(data)
  if(npc == 1L | is.data.frame(data)) data <- as.matrix(data)

  ## compute unique non-NA values
  data_unique <- as.vector(na.omit(unique(as.vector(data))))

  ## compute number of subjects and objects
  nsubj <- NROW(data)
  nobj <- 0.5 + sqrt(0.25 + if(ordered) npc else 2 * npc)
  
  ## covariates should be data.frame(-like)
  if(!is.null(covariates)) covariates <- as.data.frame(covariates)
  
  ## sanity checks
  stopifnot(is.matrix(data), nobj >= 2, isTRUE(all.equal(nobj, round(nobj))),
    isTRUE(all.equal(data_unique, round(data_unique))))
  if(!is.null(covariates)) stopifnot(nrow(covariates) == nobj)

  ## coerce to integer and set dimnames
  cnam <- which(upper.tri(diag(nobj)), arr.ind = TRUE)
  if(ordered) cnam <- rbind(cnam, cnam[,2:1])
  cnam <- apply(cnam, 1, paste, collapse = ":")
  data <- structure(as.integer(data), .Dim = c(nsubj, npc), .Dimnames = list(rownames(data), cnam))
  data_unique <- as.integer(data_unique)

  ## process labels
  if(is.null(labels)) labels <- letters[1:nobj] else {
    if(length(labels) != nobj) stop("length of labels does not match number of objects")
  }

  ## process mscale
  if(is.null(mscale)) {
    mscale <- if(length(data_unique) < 1) c(-1, 1) else {
      mscale <- max(abs(data_unique))
      if(!(mscale > 0)) mscale <- 1
      mscale <- (-mscale):mscale
      if(!(0 %in% data_unique)) mscale <- mscale[-(length(mscale)+1)/2]
      mscale
      }
  } else {
    if(!all(data_unique %in% mscale)) stop("mscale does not match data")
    mscale <- as.integer(sort(mscale))
    if(max(abs(mscale)) <= 0) stop("mscale needs to have non-zero elements")
    if(abs(head(mscale, 1)) != tail(mscale, 1)) stop("mscale must by symmetric")
  }

  ## process covariates
  if(!is.null(covariates)) rownames(covariates) <- seq_along(labels)

  rval <- data
  class(rval) <- "paircomp"
  attributes(rval) <- c(attributes(rval), list(
    labels = labels,
    mscale = mscale,
    ordered = ordered,
    covariates = covariates))
  class(rval) <- "paircomp"
  return(rval)
}



## Methods: format() handles character formatting
##   utilized in print() and as.character() method.

format.paircomp <- function(x, sep = ", ", brackets = TRUE,
  abbreviate = NULL, width = getOption("width") - 7, ...)
{
  ## process brackets
  if(is.null(brackets)) brackets <- TRUE
  if(is.logical(brackets)) brackets <- if(brackets) c("{", "}") else ""
  brackets <- rep(as.character(brackets), length.out = 2)

  ## set up comparison symbols
  mscale <- as.integer(sort(unique(c(mscale(x), 0))))
  mscale_max <- tail(mscale, 1)
  mscale_symbol <- if(mscale_max > 2L) {
    paste(ifelse(abs(mscale) > 1L, abs(mscale), ""), c("<", "=", ">")[sign(mscale) + 2L], sep = "")
  } else if(mscale_max == 2L) {
    c("<<", "<", "=", ">", ">>")[mscale + 3L]
  } else {
    c("<", "=", ">")[mscale + 2L]
  }

  ## process labels
  lab <- labels(x)
  ## abbreviate
  if(is.null(abbreviate)) abbreviate <- TRUE
  if(is.logical(abbreviate)) {
    nlab <- max(nchar(lab))
    abbreviate <- if(abbreviate) as.numeric(cut(nlab, c(-Inf, 1.5, 4.5, 7.5, Inf))) else nlab
  }
  lab <- abbreviate(lab, abbreviate)  
  ## expand
  ix <- which(upper.tri(diag(length(lab))), arr.ind = TRUE)
  lab1 <- lab[ix[,1]]
  lab2 <- lab[ix[,2]]
  if(attr(x, "ordered")) {
    lab1 <- c(lab1, lab2)
    lab2 <- c(lab2, lab1)
  }
  
  pc <- as.matrix(x)
  pc <- pc + mscale_max + 1
  pc <- apply(pc, 1, function(x) paste(
    ifelse(is.na(x), "NA", paste(lab1, mscale_symbol[x], lab2, sep = " ")),
    collapse = ", "))
  ## Handle NAs differently? Maybe via "a ? b"?

  ## check width
  if(is.null(width)) width <- TRUE
  if(is.logical(width)) width <- if(width) getOption("width") - 7 else Inf
  width <- width - sum(nchar(brackets))
  wi <- nchar(pc) > width
  if(any(wi)) {
    pc[wi] <- paste(substr(pc[wi], 1, width - 3), "...", sep = "")
  }
  
  ## add brackets
  pc <- paste(brackets[1], pc, brackets[2], sep = "")

  ## assure names (if any)
  names(pc) <- rownames(unclass(x))

  return(pc)
}

print.paircomp <- function(x, quote = FALSE, ...)
{
  print(format(x, ...), quote = quote)
  invisible(x)
}



## Basic summaries:
## str() just shows structure (no data),
## summary() shows frequency table for each comparison.

str.paircomp <- function(object, width = getOption("width") - 7, ...)
{
  rval <- sprintf(" Paired comparisons from %d subjects for %d objects: %s.\n",
    nrow(unclass(object)), length(labels(object)), paste(labels(object), collapse = ", "))
  
  if(is.null(width)) width <- TRUE
  if(is.logical(width)) width <- if(width) getOption("width") - 7 else Inf
  if(nchar(rval) > width) rval <- paste(substr(rval, 1, width - 3), "...\n", sep = "")
  
  cat(rval)
  invisible(NULL)
}

summary.paircomp <- function(object, abbreviate = FALSE, decreasing = TRUE, ...)
{
  ## data
  dat <- as.matrix(object)
  
  ## process labels
  lab <- labels(object)
  ## abbreviate
  if(is.null(abbreviate)) abbreviate <- TRUE
  if(is.logical(abbreviate)) {
    nlab <- max(nchar(lab))
    abbreviate <- if(abbreviate) as.numeric(cut(nlab, c(-Inf, 1.5, 4.5, 7.5, Inf))) else nlab
  }
  lab <- abbreviate(lab, abbreviate)  

  ## rownames
  ix <- which(upper.tri(diag(length(lab))), arr.ind = TRUE)
  lab1 <- lab[ix[,1]]
  lab2 <- lab[ix[,2]]
  if(attr(object, "ordered")) {
    lab1 <- c(lab1, lab2)
    lab2 <- c(lab2, lab1)
  }
  rnam <- paste(format(lab1), ":", format(lab2))

  ## colnames
  mscale <- mscale(object)
  mscale_max <- tail(mscale, 1)
  cnam <- if(mscale_max > 2L) {
    paste(ifelse(abs(mscale) > 1L, abs(mscale), ""), c("<", "=", ">")[sign(mscale) + 2L], sep = "")
  } else if(mscale_max == 2L) {
    c("<<", "<", "=", ">", ">>")[mscale + 3L]
  } else {
    c("<", "=", ">")[mscale + 2L]
  }
  if(decreasing) cnam <- rev(cnam)
  if(any(is.na(dat))) cnam <- c(cnam, "NA's")

  rval <- t(apply(dat, 2, function(x) table(factor(x, levels = mscale(object)))))
  if(decreasing) rval <- rval[,ncol(rval):1]
  if(any(is.na(dat))) rval <- cbind(rval, apply(dat, 2, function(x) sum(is.na(x))))
  dimnames(rval) <- list(rnam, cnam)
  
  rval
}



## extract and combine observations

length.paircomp <- function(x) nrow(unclass(x))

"[.paircomp" <- function(x, i, ...) {
  xattr <- attributes(x)[c("labels", "mscale", "ordered", "covariates")]
  x <- unclass(x)[i,,drop=FALSE]
  attributes(x) <- c(attributes(x), xattr)
  structure(x, class = "paircomp")
}

c.paircomp <- function(...)
{
  args <- list(...)

  check_list <- function(x) {
    if(length(x) < 2) return(TRUE)
    x1 <- x[[1]]
    all(sapply(2:length(x), function(i) identical(x[[i]], x1)))
  }
  if(!check_list(lapply(args, attr, "labels"))) stop("objects have different labels")
  if(!check_list(lapply(args, attr, "mscale"))) stop("objects have different mscales")
  if(!check_list(lapply(args, attr, "ordered"))) stop("objects have differently ordered")
  if(!check_list(lapply(args, attr, "covariates"))) stop("objects have different covariates")

  rval <- args[[1]]
  xattr <- attributes(rval)[c("labels", "mscale", "ordered", "covariates")]
  rval <- do.call("rbind", lapply(args, function(x) as.matrix(x)))
  attributes(rval) <- c(attributes(rval), xattr)
  structure(rval, class = "paircomp")
}

rep.paircomp <- function(x, ...) {
  ix <- if(length(x) > 0) 1:length(x) else 0
  x[rep(ix, ...)]
}


## mscale(): new generic with extractor method for paircomp

mscale <- function(object, ...) UseMethod("mscale")

mscale.paircomp <- function(object, ...) attr(object, "mscale")



## covariates(): new generic with extractor method for paircomp

covariates <- function(object, ...) UseMethod("covariates")

"covariates<-" <- function(object, value) UseMethod("covariates<-")

covariates.paircomp <- function(object, ...) {
  if(is.null(attr(object, "covariates"))) return(NULL)
  rval <- attr(object, "covariates")
  rownames(rval) <- labels(object)
  return(rval)
}

"covariates<-.paircomp" <- function(object, value) {
  if(!is.data.frame(value)) {
    dval <- as.data.frame(value)
    if(ncol(dval) == 1) names(dval) <- paste(deparse(substitute(value), width.cutoff = 500), collapse = " ")
    value <- dval
  }
  if(nrow(value) != length(labels(object))) stop("Number of rows in covariates does not match number of objects.")
  rownames(value) <- 1:length(labels(object))
  attr(object, "covariates") <- value
  return(object)
}


## names() queries/sets subject names,
## labels() queries/sets object names, (new generic labels<-)
## levels() is an alias for labels() (but discouraged to avoid confusion with mscale)

names.paircomp <- function(x) rownames(unclass(x))

"names<-.paircomp" <- function(x, value) {
  x <- unclass(x)
  rownames(x) <- value
  structure(x, class = "paircomp")
}

"labels<-" <- function(object, value) UseMethod("labels<-")
  
labels.paircomp <- function(object, ...) attr(object, "labels")

"labels<-.paircomp" <- function(object, value) {
  if(!(length(value) == length(attr(object, "labels")))) stop("length of labels does not match")
  attr(object, "labels") <- value
  object
}

levels.paircomp <- function(x, ...) attr(object, "labels")

"levels<-.paircomp" <- function(x, value) {
  labels(x) <- value
  x
}



## reorder() method enables selection of subsets
## of objects to be compared and/or re-ordering of
## objects.
## subset() currently just calls reorder(), in addition
## selection of subsets of subjects would be nice to have.

reorder.paircomp <- function(x, labels, ...)
{
  xlab <- labels(x)
  if(missing(labels)) labels <- xlab
  if(is.character(labels)) labels <- sapply(labels, match.arg, choices = xlab)
  if(is.character(labels)) labels <- match(labels, xlab)
  if(length(labels) < 2) stop("Need to compare at least two objects.")
  
  ## select relevant columns
  nlab <- labels
  ix <- which(upper.tri(diag(length(xlab))), arr.ind = TRUE)
  wi <- apply(ix, 1, function(z) all(z %in% nlab))
  ndat <- if(attr(x, "ordered")) as.matrix(x)[, c(wi, wi), drop = FALSE] else as.matrix(x)[, wi, drop = FALSE]
  
  ## re-order comparisons (if necessary)
  if(!identical(nlab, sort(nlab))) {
    ## set up indexes
    ix <- ix[wi,, drop = FALSE]
    nix <- which(upper.tri(diag(length(nlab))), arr.ind = TRUE)
    nix <- matrix(nlab[nix], ncol = 2)
    ## match
    wi <- apply(nix, 1, function(x) x[1] > x[2])
    nix[wi,] <- nix[wi, 2:1]
    onam <- apply(ix, 1, paste, collapse = ":")
    nnam <- apply(nix, 1, paste, collapse = ":")
    ord <- match(nnam, onam)
    
    ## re-order
    if(attr(x, "ordered")) {
      ndat <- ndat[, c(ord, ord), drop = FALSE]
      ndat[, c(wi, wi)] <- -1 * ndat[, c(wi, wi)]
    } else {
      ndat <- ndat[, ord, drop = FALSE]
      ndat[, wi] <- -1 * ndat[, wi]
    }    
  }

  ## call constructor to setup proper new paircomp object
  paircomp(ndat,
    labels = attr(x, "labels")[nlab], mscale = mscale(x),
    ordered = attr(x, "ordered"), covariates = attr(x, "covariates")[nlab,,drop=FALSE])
}

subset.paircomp <- function(x, subset, select, ...)
{
  if(!missing(subset)) stop("'subset' argument not yet implemented.")
  ## FIXME: This should enable the user to select all observations with, say:
  ## "a > b" or "a : b > -1" etc.

  reorder(x, labels = select, ...)
}


## coercion functions, all fairly straightforward

as.character.paircomp <- function(x, ...) {
  format(x, abbreviate = FALSE, width = FALSE, ...)
}

as.integer.paircomp <-
as.matrix.paircomp <- function(x, ...) {
  rval <- unclass(x)
  attr(rval, "labels") <- attr(rval, "mscale") <- attr(rval, "ordered") <- attr(rval, "covariates") <- NULL

  ## colnames
  lab <- attr(x, "labels")
  ix <- which(upper.tri(diag(length(lab))), arr.ind = TRUE)
  lab1 <- lab[ix[,1]]
  lab2 <- lab[ix[,2]]
  if(attr(x, "ordered")) {
    lab1 <- c(lab1, lab2)
    lab2 <- c(lab2, lab1)
  }
  colnames(rval) <- paste(lab1, ":", lab2, sep = "")
  return(rval)  
}

as.double.paircomp <- function(x, ...) {
  rval <- as.matrix(x, ...)
  rval[] <- as.double(rval)
  return(rval)
}

as.data.frame.paircomp <- function(x, ...) {
  rval <- as.data.frame(numeric(length(x)), ...)
  rval[[1]] <- x
  names(rval)[1] <- paste(deparse(substitute(x), width.cutoff = 500), collapse = " ")
  return(rval)
}


## FIXME: na.* methods needed

