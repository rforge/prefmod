\name{btReg}
\alias{btReg}
\alias{btReg.fit}

\alias{print.btReg}
\alias{reweight.btReg}
\alias{summary.btReg}
\alias{print.summary.btReg}

\alias{coef.btReg}
\alias{worth.btReg}
\alias{deviance.btReg}
\alias{logLik.btReg}
\alias{vcov.btReg}
\alias{bread.btReg}
\alias{estfun.btReg}
\alias{plot.btReg}


\title{Model-Generating Function for Bradley-Terry Models}

\description{
  \code{btReg} is a wrapper function that creates a \code{"StatModel"} object
  with certain fitting parameters passed on to \code{btReg.fit} for fitting
  simple Bradley-Terry models.
}

\usage{
btReg(type = "loglin", ref = NULL, undecided = NULL, position = NULL)

btReg.fit(y, weights = NULL, type = c("loglin", "logit"), ref = NULL,
  undecided = NULL, position = NULL, \dots)
}

\arguments{
  \item{y}{paircomp object with the response.}
  \item{weights}{an optional vector of weights, interpreted as case weights (integer only).}
  \item{type}{character. Should an auxiliary log-linear Poisson model or logistic
    binomial be employed for estimation? The latter is only available if not undecided
    effects are estimated.}
  \item{ref}{character or numeric. Which object parameter should be the reference
    category, i.e., constrained to zero?}
  \item{undecided}{logical. Should an undecided parameter be estimated?}
  \item{position}{logical. Should a position effect be estimated?}
  \item{\dots}{further arguments passed to functions.}
}

\details{
  \code{btReg} is the S4 interface for \code{btReg.fit} so that it
  can be used in \code{\link[party]{mob}}. The user does not have to
  call this directly but can simply use the \code{\link{bttree}}
  interface.
}

\value{
  An object of class \code{"StatModel"} that fits Bradley-Terry models
  with the specified arguments.
}

\seealso{\code{\link{bttree}}}

\keyword{tree}
