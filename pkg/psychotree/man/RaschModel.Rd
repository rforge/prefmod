\name{RaschModel}

\alias{RaschModel}
\alias{reweight.RaschModel}
\alias{bread.RaschModel}
\alias{estfun.RaschModel}


\title{Model-Generating Function for Rasch Models}

\description{
  \code{RaschModel} is a wrapper function that creates a \code{"StatModel"} object
  with certain fitting parameters passed on to \code{\link[psychotools]{RaschModel.fit}} for fitting
  simple Rasch models.
}

\usage{
RaschModel(gradtol = 1e-6, deriv = c("sum", "diff", "numeric"),
  hessian = TRUE, iterlim = 100L)
}

\arguments{
  \item{deriv}{character. Which type of derivatives should be used for computing
    gradient and Hessian matrix? Analytical with sum algorithm (\code{"sum"}),
    analytical with difference algorithm (\code{"diff"}, faster but numerically unstable),
    or numerical.}
  \item{hessian}{logical. Should the Hessian of the final model be computed?
    If set to \code{FALSE}, the \code{vcov} method can only return \code{NA}s
    and consequently no standard errors or tests are available in the \code{summary}.}
  \item{gradtol, iterlim, \dots}{further arguments passed to \code{\link[stats]{nlm}}.}
}

\details{
  \code{RaschModel} is the S4 interface for \code{\link[psychotools]{RaschModel.fit}} so that it
  can be used in \code{\link[party]{mob}}. The user does not have to
  call this directly but can simply use the \code{\link{raschtree}}
  interface.
}

\value{
  \code{RaschModel} returns an S4 object of class \code{"StatModel"} that fits Rasch models
  with the specified arguments.
}

\seealso{\code{\link{raschtree}}, \code{\link[psychotools]{RaschModel.fit}}}

\keyword{tree}
