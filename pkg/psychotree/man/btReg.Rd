\name{btReg}
\alias{btReg}
\alias{print.btReg}
\alias{reweight.btReg}

\title{Model-Generating Function for Bradley-Terry Models}

\description{
  A wrapper function that creates a \code{"StatModel"} object
  with certain fitting parameters passed on to \code{\link[prefmod2]{btreg.fit}}.
}

\usage{
btReg(type = "loglin", ref = NULL, undecided = NULL, position = NULL)
}

\arguments{
  \item{type, ref, undecided, position}{arguments for the Bradley-Terry
    model passed on to \code{\link[prefmod2]{btreg}}.}
}

\details{
  \code{btReg} is the S4 interface for \code{btreg.fit} so that it
  can be used in \code{\link[party]{mob}}. The user does not have to
  call this directly but can simply use the \code{\link{bttree}}
  interface.
}

\value{
  An object of class \code{"StatModel"} that fits Bradley-Terry models
  with the specified arguments.
}

\seealso{\code{\link{bttree}}, \code{\link[prefmod2]{btreg.fit}}}

\keyword{tree}
