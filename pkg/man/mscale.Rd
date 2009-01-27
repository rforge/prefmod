\name{mscale}

\alias{mscale}

\title{Extract Measurement Scale}

\description{
  A generic function for extracting the measurement scale from an object.
}

\usage{
  mscale(object, \dots)
}

\arguments{
  \item{object}{an object.}
  \item{\dots}{arguments passed to methods.}
}

\examples{
## method for "paircomp" data
pc <- paircomp(rbind(
  c(1,  1,  1), # a > b, a > c, b > c
  c(1,  1, -1), # a > b, a > c, b < c
  c(1, -1, -1), # a > b, a < c, b < c
  c(1,  1,  1)))
mscale(pc)
}

\keyword{classes}
