\name{PCModel}
\alias{PCModel}
\alias{reweight.PCModel}
\alias{estfun.PCModel}

\title{Model-Generating Function for Partial Credit Models}

\description{
  \code{PCModel} is a wrapper function that creates a \code{"StatModel"} object
  with certain fitting parameters passed on to \code{\link[psychotools]{PCModel.fit}} for fitting
  partial credit models.
}

\usage{
PCModel(nullcats = c("keep", "downcode", "ignore"), reltol = 1e-10,
  deriv = c("sum", "diff"), hessian = TRUE, maxit = 100L)
}

\arguments{
  \item{nullcats}{character, specifying how items with null categories,
    i.e., categories not used, are treated (see details in \code{\link[psychotools]{PCModel.fit}}).}

  \item{deriv}{character. If "sum" (the default), the first derivatives
    of the elementary symmetric functions are calculated with the sum
    algorithm. Otherwise ("diff") the difference algorithm (faster but
    numerically unstable) is used.}

  \item{hessian}{logical. Should the Hessian of the final model be computed?
    If set to \code{FALSE}, the \code{vcov} method can only return \code{NA}s
    and consequently no standard errors or tests are available in the \code{summary}.}

  \item{reltol, maxit, \dots}{further arguments passed to \code{\link[stats]{optim}}.}
}

\details{
  \code{PCModel} is the S4 interface for \code{\link[psychotools]{PCModel.fit}} so that it
  can be used in \code{\link[party]{mob}}. The user does not have to
  call this directly but can simply use the \code{\link{pctree}}
  interface.
}

\value{
  \code{PCModel} returns an S4 object of class \code{"StatModel"} that
  fits partial credit models with the specified arguments.
}

\seealso{\code{\link{pctree}}, \code{\link[psychotools]{PCModel.fit}}}

\keyword{tree}
