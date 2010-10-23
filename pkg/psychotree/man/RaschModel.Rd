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


\title{Model-Generating Function for Rasch Models}

\description{
  \code{RaschModel} is a wrapper function that creates a \code{"StatModel"} object
  with certain fitting parameters passed on to \code{RaschModel.fit} for fitting
  simple Rasch models.
}

\usage{
RaschModel(gradtol = 1e-6, deriv = c("sum", "diff", "numeric"))

RaschModel.fit(y, weights = NULL, start = NULL, gradtol = 1e-6,
  deriv = c("sum", "diff", "numeric"), \dots)
}

\arguments{
  \item{y}{object that can be coerced (via \code{\link[base]{as.matrix}}) to binary 0/1 matrix.}
  \item{weights}{an optional vector of weights, interpreted as case weights (integer only).}
  \item{start}{an optional vector of starting values.}
  \item{deriv}{character. Which type of derivatives should be used for computing
    gradient and Hessian matrix? Analytical with sum algorithm (\code{"sum"}),
    analytical with difference algorithm (\code{"diff"}, faster but numerically unstable),
    or numerical.}
  \item{gradtol, \dots}{further arguments passed to \code{\link[stats]{nlm}}.}
}

\details{
  \code{RaschModel} is the S4 interface for \code{RaschModel.fit} so that it
  can be used in \code{\link[party]{mob}}. The user does not have to
  call this directly but can simply use the \code{\link{raschtree}}
  interface.
}

\value{
  \code{RaschModel} returns an S4 object of class \code{"StatModel"} that fits Rasch models
  with the specified arguments. \code{RaschModel.fit} returns an S3 object of class
  \code{"RaschModel"}, i.e., a list with components as follows.
  \item{coefficients}{estimated item difficulty parameters (without first item
    parameter which is always constrained to be 0),}
  \item{vcov}{covariance matrix of the parameters in the model,}
  \item{loglik}{log-likelihood of the fitted model,}
  \item{df}{number of estimated parameters,}
  \item{data}{the original data supplied (excluding columns without variance),}
  \item{weights}{the weights used (if any),}
  \item{n}{number of observations (with non-zero weights),}
  \item{items}{status indicator (0, 0/1, 1) of all original items,}
  \item{na}{logical indicating whether the data contains NAs,}
  \item{elementary_symmetric_functions}{List of elementary symmetric functions
    for estimated parameters (up to order 2; or 1 in case of numeric derivatives),}
  \item{nlm_code}{convergence code from \code{nlm},}
  \item{iterations}{number of iterations used by \code{nlm},}
  \item{gradtol}{tolerance passed to \code{nlm}.}
}

\seealso{\code{\link{raschtree}}}

\keyword{tree}
