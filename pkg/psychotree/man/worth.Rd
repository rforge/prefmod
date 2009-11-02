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

\seealso{\code{\link{bttree}}}

\examples{
data("CEMSChoice", package = "psychotree")
cems_tree <- bttree(preference ~ french + spanish + italian + study + work + gender + intdegree,
  data = CEMSChoice, minsplit = 5, ref = "London")
worth(cems_tree)
}

\keyword{classes}
