\encoding{UTF-8}
\name{summary.pattNPML}
\alias{summary.pattNPML}
\alias{print.pattNPML}
\alias{BIC.pattNPML}
\title{Summarizing finite mixture regression fits}
\description{%
  These functions are the \code{summary}, \code{print}, and \code{BIC} methods for objects of type \code{pattNPML}.
}
\usage{
\method{summary}{pattNPML}(object, digits = max(3, getOption("digits") - 3), \dots)
\method{print}{pattNPML}(x, digits = max(3, getOption("digits") - 3), \dots)
}
\arguments{%
  \item{object}{a fitted object of class \code{pattNPML}.}
  \item{x}{a fitted object of class \code{pattNPML}.}
  \item{digits}{number of digits; applied on various displayed quantities.}
  \item{\dots}{further arguments, which will mostly be ignored.}
}
\note{%
  The \code{summary} and \code{print} methods are adapted versions from the \pkg{npmlreg} package.
}
\seealso{%
  \code{\link{pattnpml.fit}}%
}
\keyword{models}
\keyword{regression}
