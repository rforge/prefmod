\name{print.pattMod}
\alias{print.pattMod}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Print methods for pattern models}
\description{Print method for objects of class \code{pattMod}.}
\usage{
\method{print}{pattMod}(x, ...)
}
\arguments{
  \item{x}{Object of class \code{pattMod}.}
  \item{...}{Further arguments to be passed to or from other methods. They are ignored in this function.}
}
\details{
  This print method generates output for fitted pattern models, i.e. for models.
  of class \code{pattMod}. The functions \code{\link{pattPC.fit}},
  \code{\link{pattR.fit}}, and \code{\link{pattPC.fit}} produce such objects.
}
\author{Reinhold Hatzinger}
\examples{
data(salad)
res<-pattR.fit(salad, nitems=4)
print(res)
}
\keyword{models}
