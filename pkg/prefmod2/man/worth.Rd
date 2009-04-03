\name{worth}

\alias{worth}

\title{Extract Worth Parameters}

\description{
  Generic functions for extracting worth parameters from paired comparison models.
}

\usage{
  worth(object, \dots)
}

\arguments{
  \item{object}{an object.}
  \item{\dots}{arguments passed to methods.}
}

\seealso{\code{\link{btreg}}}

\examples{
data("CEMSChoice", package = "prefmod2")
cems <- btreg(preference ~ 1, data = CEMSChoice)
worth(cems)
}

\keyword{classes}
