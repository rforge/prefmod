\name{RaschModel}
\alias{RaschModel}
\alias{RaschModel.fit}

\alias{print.RaschModel}
\alias{reweight.RaschModel}
\alias{summary.RaschModel}
\alias{print.summary.RaschModel}

\alias{coef.RaschModel}
\alias{worth.RaschModel}
\alias{logLik.RaschModel}
\alias{vcov.RaschModel}
\alias{bread.RaschModel}
\alias{estfun.RaschModel}
\alias{plot.RaschModel}


\title{Model-Generating Function for Rasch Models}

\description{
  \code{RaschModel} is a wrapper function that creates a \code{"StatModel"} object
  with certain fitting parameters passed on to \code{RaschModel.fit} for fitting
  simple Rasch models.
}

\usage{
RaschModel(gradtol = 1e-12, hessian = FALSE)

RaschModel.fit(y, weights = NULL, start = NULL,
  gradtol = 1e-12, check.analyticals = FALSE, hessian = FALSE, \dots)
}

\arguments{
  \item{y}{object that can be coerced (via \code{\link[base]{as.matrix}}) to binary 0/1 matrix.}
  \item{weights}{an optional vector of weights, interpreted as case weights (integer only).}
  \item{start}{an optional vector of starting values.}
  \item{hessian}{logical. Should the numerical Hessian be used for inference? (Default
    is to use analytical Hessian.)}
  \item{gradtol, check.analyticals, \dots}{further arguments passed to \code{\link[stats]{nlm}}.}
}

\details{
  \code{RaschModel} is the S4 interface for \code{RaschModel.fit} so that it
  can be used in \code{\link[party]{mob}}. The user does not have to
  call this directly but can simply use the \code{\link{raschtree}}
  interface.
}

\value{
  An object of class \code{"StatModel"} that fits Rasch models
  with the specified arguments.
}

\seealso{\code{\link{raschtree}}}

\keyword{tree}
