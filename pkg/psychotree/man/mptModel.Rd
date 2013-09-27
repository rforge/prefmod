\name{mptModel}

\alias{mptModel}
\alias{reweight.mptModel}
\alias{bread.mptModel}
\alias{estfun.mptModel}


\title{Model-Generating Function for Multinomial Processing Tree Models}

\description{
  \code{mptModel} is a wrapper function that creates a \code{"StatModel"}
  object with certain fitting parameters passed on to
  \code{\link[psychotools]{mptModel.fit}} for fitting simple MPT models.
}

\usage{
mptModel(mptform = NULL, mptstruc = NULL, maxit = 1000)
}

\arguments{
  \item{mptform}{formula representation of the MPT model structure, usually
    returned by \code{\link{mpttree}}.}
  \item{mptstruc}{list of structural constants, usually returned by
    \code{\link{mpttree}}.}
  \item{maxit}{maximum number of iterations, passed to
    \code{\link[psychotools]{mptModel.fit}}.}
}

\details{
  \code{mptModel} is the S4 interface for
  \code{\link[psychotools]{mptModel.fit}} so that it can be used in
  \code{\link[party]{mob}}. The user does not have to call this directly but
  can simply use the \code{\link{mpttree}} interface.
}

\value{
  \code{mptModel} returns an S4 object of class \code{"StatModel"} that fits
  multinomial processing tree (MPT) models with the specified arguments.
}

\seealso{\code{\link{mpttree}}, \code{\link[psychotools]{mptModel.fit}}}

\keyword{tree}
