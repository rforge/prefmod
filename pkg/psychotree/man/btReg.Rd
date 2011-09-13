\name{btReg}

\alias{btReg}
\alias{reweight.btReg}
\alias{bread.btReg}
\alias{estfun.btReg}


\title{Model-Generating Function for Bradley-Terry Models}

\description{
  \code{btReg} is a wrapper function that creates a \code{"StatModel"} object
  with certain fitting parameters passed on to \code{\link[psychotools]{btReg.fit}} for fitting
  simple Bradley-Terry models.
}

\usage{
btReg(type = "loglin", ref = NULL, undecided = NULL, position = NULL)
}

\arguments{
  \item{type}{character. Should an auxiliary log-linear Poisson model or logistic
    binomial be employed for estimation? The latter is only available if not undecided
    effects are estimated.}
  \item{ref}{character or numeric. Which object parameter should be the reference
    category, i.e., constrained to zero?}
  \item{undecided}{logical. Should an undecided parameter be estimated?}
  \item{position}{logical. Should a position effect be estimated?}
}

\details{
  \code{btReg} is the S4 interface for \code{\link[psychotools]{btReg.fit}} so that it
  can be used in \code{\link[party]{mob}}. The user does not have to
  call this directly but can simply use the \code{\link{bttree}}
  interface.
}

\value{
  An object of class \code{"StatModel"} that fits Bradley-Terry models
  with the specified arguments.
}

\seealso{\code{\link{bttree}}, \code{\link[psychotools]{btReg.fit}}}

\keyword{tree}
