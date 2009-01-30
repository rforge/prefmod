\name{print.pattW}
\alias{print.pattW}
\title{Print worth matrix for pattern models}
\description{Print method for objects of class \code{pattW}.}
\usage{
\method{print}{pattW}(x, ...)
}
\arguments{
  \item{x}{Object of class \code{pattW}.}
  \item{...}{Further arguments to be passed to or from other methods. They are ignored in this function.}
}
\details{
  This print method prints the output from \code{link{patt.worth}}.
}
\author{Reinhold Hatzinger}
\examples{
data(salad)
res<-pattR.fit(salad, nitems=4)
wres<-patt.worth(res)
print(wres)
}
\keyword{models}
