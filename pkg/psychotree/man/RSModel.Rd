\name{RSModel}

\alias{RSModel}
\alias{reweight.RSModel}
\alias{estfun.RSModel}


\title{Model-Generating Function for Rating Scale Models}


\description{

  \code{RSModel} is a wrapper function that creates a \code{"StatModel"} object
  with certain fitting parameters passed on to \code{\link[psychotools]{RSModel.fit}} for fitting
  rating scale models.

}


\usage{

RSModel(reltol = 1e-10, deriv = c("sum", "diff"),
  hessian = TRUE, maxit = 100L)

}

\arguments{

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

  \code{RSModel} is the S4 interface for \code{\link[psychotools]{RSModel.fit}} so that it
  can be used in \code{\link[party]{mob}}. The user does not have to
  call this directly but can simply use the \code{\link{rstree}}
  interface.

}


\value{

  \code{RSModel} returns an S4 object of class \code{"StatModel"} that
  fits rating scale models with the specified arguments.

}


\seealso{\code{\link{rstree}}, \code{\link[psychotools]{RSModel.fit}}}


\keyword{tree}
