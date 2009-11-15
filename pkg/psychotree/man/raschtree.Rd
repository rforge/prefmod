\name{raschtree}

\alias{raschtree}

\alias{plot.raschtree}
\alias{print.raschtree}
\alias{summary.raschtree}
\alias{sctest.raschtree}

\alias{deviance.raschtree}
\alias{logLik.raschtree}
\alias{coef.raschtree}
\alias{worth.raschtree}
\alias{weights.raschtree}

\title{Rasch Tree Models}

\description{
  Recursive partitioning based on Rasch models.
}

\usage{
raschtree(formula, data, na.action = na.pass, minsplit = 10, \dots)
}

\arguments{
  \item{formula}{A symbolic description of the model to be fit. This
    should be of type \code{y ~ x1 + x2}
    where \code{y} should be a binary 0/1 matrix
    and \code{x1} and \code{x2} are used as partitioning variables.}
  \item{data}{an optional data frame containing the variables in the model.}
  \item{na.action}{A function which indicates what should happen when the data
    contain \code{NA}s, defaulting to \code{\link[stats]{na.pass}}.}
  \item{minsplit, \dots}{arguments passed to \code{\link[party]{mob_control}}.}
}

\details{
  Rasch tree models are an application of model-based recursive partitioning
  (implemented in \code{\link[party]{mob}}) to Rasch models
  (implemented in \code{\link{RaschModel}}).
  
  Various methods are provided for \code{"raschtree"} objects, most of them
  inherit their behavior from \code{"mob"} objects (e.g., \code{print}, \code{summary},
  etc.). \code{worth} behaves analogously to \code{coef} and extracts the
  item parameters from the Rasch models in the nodes of the tree. The \code{plot}
  method employs the \code{\link{node_raschplot}} panel-generating function.
}

\value{
  An object of S3 class \code{"raschtree"} which is a list containing only
  a single element of S4 class \code{"mob"} (because this is currently not
  exported).
}

\seealso{\code{\link[party]{mob}}, \code{\link{RaschModel}}, \code{\link{bttree}}}

\keyword{tree}
