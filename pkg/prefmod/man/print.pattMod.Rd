\encoding{UTF-8}
\name{print.pattMod}
\alias{print.pattMod}
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
  of class \code{pattMod}. The functions \code{\link{pattPC.fit}}, \code{\link{pattR.fit}},
  \code{\link{pattL.fit}}, and  \code{\link{pattLrep.fit}} produce such objects.
}
\author{Reinhold Hatzinger}
\examples{
res <- pattR.fit(salad, nitems = 4)
res
}
\keyword{models}
