paircomp <- function(data, labels = NULL, scale = NULL, ordered = FALSE)
{
  ## data should be matrix(-like)
  npc <- NCOL(data)
  if(npc == 1L | is.data.frame(data)) data <- as.matrix(data)

  ## compute unique non-NA values
  data_unique <- as.vector(na.omit(unique(as.vector(data))))

  ## compute number of subjects and objects
  nsubj <- NROW(data)
  nobj <- 0.5 + sqrt(0.25 + if(ordered) npc else 2 * npc)
  
  ## sanity checks
  stopifnot(is.matrix(data), nobj >= 2, isTRUE(all.equal(nobj, round(nobj))),
    isTRUE(all.equal(data_unique, round(data_unique))))

  ## coerce to integer
  data <- structure(as.integer(data), .Dim = c(nsubj, npc), .Dimnames = list(rownames(data), NULL))
  data_unique <- as.integer(data_unique)

  ## process labels
  if(is.null(labels)) labels <- letters[1:nobj] else {
    if(length(labels) != nobj) stop("length of labels does not match number of objects")
  }

  ## process scale
  if(is.null(scale)) {
    scale <- if(length(data_unique) < 1) c(-1, 1) else {
      scale <- max(abs(data_unique))
      if(!(scale > 0)) scale <- 1
      scale <- (-scale):scale
      if(!(0 %in% data_unique)) scale <- scale[-(length(scale)+1)/2]
      scale
      }
  } else {
    if(!all(data_unique %in% scale)) stop("scale does not match data")
    scale <- as.integer(sort(scale))
  }

  rval <- list(
    data = data,
    labels = labels,
    scale = scale,
    ordered = ordered)
  class(rval) <- "paircomp"
  return(rval)
}

## FIXME: scale?

format.paircomp <- function(x, sep = ", ", brackets = TRUE,
  abbreviate = NULL, width = getOption("width") - 7, ...)
{
  ## FIXME: handle NAs!
  ## write: "NA" or "a ? b"

  ## process brackets
  if(is.null(brackets)) brackets <- TRUE
  if(is.logical(brackets)) brackets <- if(brackets) c("{", "}") else ""
  brackets <- rep(as.character(brackets), length.out = 2)

  ## set up comparison symbols
  scale <- as.integer(sort(unique(c(x$scale, 0))))
  scale_max <- tail(scale, 1)
  scale_symbol <- if(scale_max > 2L) {
    paste(ifelse(abs(scale) > 1L, abs(scale), ""), c("<", "=", ">")[sign(scale) + 2L], sep = "")
  } else if(scale_max == 2L) {
    c("<<", "<", "=", ">", ">>")[scale + 3L]
  } else {
    c("<", "=", ">")[scale + 2L]
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
  if(x$ordered) {
    lab1 <- c(lab1, lab2)
    lab2 <- c(lab2, lab1)
  }
  
  pc <- x$data
  pc <- pc + scale_max + 1
  pc <- apply(pc, 1, function(x) paste(lab1, scale_symbol[x], lab2, sep = " ", collapse = ", "))

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
  names(pc) <- rownames(x$data)

  return(pc)
}

print.paircomp <- function(x, quote = FALSE, ...)
{
  print(format(x, ...), quote = quote)
  invisible(x)
}

"[.paircomp" <- function(x, i) {
  x$data <- x$data[i,,drop=FALSE]
  return(x)
}

c.paircomp <- function(...)
{
  args <- list(...)

  check_list <- function(x) {
    if(length(x) < 2) return(TRUE)
    x1 <- x[[1]]
    all(sapply(2:length(x), function(i) identical(x[[i]], x1)))
  }
  if(!check_list(lapply(args, function(x) x$labels))) stop("objects have different labels")
  if(!check_list(lapply(args, function(x) x$scale))) stop("objects have different scales")
  if(!check_list(lapply(args, function(x) x$ordered))) stop("objects have differently ordered")

  rval <- args[[1]]
  rval$data <- do.call("rbind", lapply(args, function(x) x$data))
  return(rval)  
}

length.paircomp <- function(x) nrow(x$data)

names.paircomp <- function(x) rownames(x$data)
"names<-.paircomp" <- function(x, value) {
  rownames(x$data) <- value
  x
}

"labels<-" <- function(object, value) UseMethod("labels<-")
  
labels.paircomp <- function(object, ...) object$labels

"labels<-.paircomp" <- function(object, value) {
  if(!(length(value) == length(object$labels))) stop("length of labels does not match")
  object$labels <- value
  object
}

levels.paircomp <- function(x, ...) object$labels

"levels<-.paircomp" <- function(x, value) {
  labels(x) <- value
  x
}

## TODO: str, summary, subset, na.*, reorder

pc <- paircomp(rbind(
  c(1,  1,  1), # a > b > c
  c(1,  1, -1), # a > c > b
  c(1, -1, -1), # c > a > b
  c(1,  1,  1)))

pc2 <- paircomp(rbind(
  c(2,  1,  0),
  c(1,  1, -1),
  c(1, -2, -1),
  c(0,  0,  0)),
  labels = c("Nordrhein-Westfalen", "Schleswig-Holstein", "Baden-Wuerttemberg"))
  
dat <- data.frame(
  x = rnorm(4),
  y = factor(c(1, 2, 1, 1), labels = c("hansi", "beppi")))
dat$pc <- pc

