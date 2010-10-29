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
raschtree(formula, data, minsplit = 10, gradtol = 1e-6,
  deriv = c("sum", "diff", "numeric"), \dots)
}

\arguments{
  \item{formula}{A symbolic description of the model to be fit. This
    should be of type \code{y ~ x1 + x2}
    where \code{y} should be a binary 0/1 matrix
    and \code{x1} and \code{x2} are used as partitioning variables.}
  \item{data}{a data frame containing the variables in the model.}
  \item{deriv}{character. Which type of derivatives should be used for computing
    gradient and Hessian matrix? Analytical with sum algorithm (\code{"sum"}),
    analytical with difference algorithm (\code{"diff"}, faster but numerically unstable),
    or numerical. Passed to \code{\link{RaschModel}}.}
  \item{gradtol}{numeric tolerance passed to \code{\link{RaschModel}}
    and on to \code{\link[stats]{nlm}}.}
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

  Rasch tree models are introduced in Strobl et al. (2010) whose analysis
  for the \code{\link{SPISA}} data is replicated in
  \code{vignette("raschtree", package = "psychotree")}. Their illustration
  employing artificial data is replicated below.
}

\references{
  Strobl, C., Kopf, J., and Zeileis, A. (2010). A New Method for Detecting
  Differential Item Functioning in the Rasch Model. Unpublished Manuscript.
}

\value{
  An object of S3 class \code{"raschtree"} which is a list containing only
  a single element of S4 class \code{"mob"} (because this is currently not
  exported from the party package).
}

\seealso{\code{\link[party]{mob}}, \code{\link{RaschModel}}, \code{\link{bttree}}}

\examples{
## artificial data
data("DIFSim", package = "psychotree")

## fit Rasch tree model
rt <- raschtree(resp ~ age + gender + motivation, data = DIFSim)
plot(rt)

## extract item parameters
coef(rt)
worth(rt)

## inspect parameter stability tests in all splitting nodes
sctest(rt, node = 1)
sctest(rt, node = 2)

## highlight items 3 and 14 with DIF
ix <- rep(1, 20)
ix[c(3, 14)] <- 2
plot(rt, ylines = 2.5,  cex = c(0.4, 0.8)[ix], 
  pch = c(19, 19)[ix], col = gray(c(0.5, 0))[ix])
}

\keyword{tree}
