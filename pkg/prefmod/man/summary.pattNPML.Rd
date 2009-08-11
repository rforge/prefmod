\name{summary.pattNPML}
\alias{summary.pattNPML}
\alias{print.pattNPML}
\title{Summarizing finite mixture regression fits}
\description{
  These functions are the \code{summary} and \code{print} methods for objects of  type
  \code{pattNPML}.
}
\usage{
\method{summary}{pattNPML}(object, digits = max(3, getOption("digits") - 3), ...)
\method{print}{pattNPML}(x, digits=max(3,getOption('digits')-3), ...)
}
\arguments{
  \item{object}{a fitted object of class \code{pattNPML}.   }
  \item{x}{a fitted object of class \code{pattNPML}.   }
  \item{digits}{number of digits; applied on various displayed quantities.}
  \item{\dots}{further arguments, which will mostly be ignored.}
}
\note{
These methods are adapted versions from the npmlreg package.
}
\seealso{ \code{\link{pattnpml.fit}}}
\keyword{ models }
\keyword{ regression }
