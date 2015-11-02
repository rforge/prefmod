\name{raschtree}
\alias{raschtree}

\alias{print.raschtree}
\alias{plot.raschtree}
\alias{predict.raschtree}
\alias{itempar.raschtree}

\title{Rasch Trees}

\description{
  Recursive partitioning (also known as trees) based on Rasch models.
}

\usage{
raschtree(formula, data, na.action,
  reltol = 1e-10, deriv = c("sum", "diff", "numeric"), maxit = 100L,
  \dots)

\method{predict}{raschtree}(object, newdata = NULL,
  type = c("probability", "cumprobability", "mode", "median", "mean",
    "category-information", "item-information", "test-information", "node"),
  personpar = 0, \dots)

\method{plot}{raschtree}(x, type = c("profile", "regions"), terminal_panel = NULL,
  tp_args = list(...), tnex = 2L, drop_terminal = TRUE, ...)
}

\arguments{
  \item{formula}{A symbolic description of the model to be fit. This
    should be of type \code{y ~ x1 + x2}
    where \code{y} should be a binary 0/1 item response matrix
    and \code{x1} and \code{x2} are used as partitioning variables.}
  \item{data}{a data frame containing the variables in the model.}
  \item{na.action}{a function which indicates what should happen when the data
    contain missing values (\code{NA}s).}
  \item{deriv}{character. Which type of derivatives should be used for computing
    gradient and Hessian matrix? Analytical with sum algorithm (\code{"sum"}),
    analytical with difference algorithm (\code{"diff"}, faster but numerically unstable),
    or numerical. Passed to \code{\link[psychotools]{raschmodel}}.}
  \item{reltol, maxit}{arguments passed via \code{\link[psychotools]{raschmodel}}
    to \code{\link[stats]{optim}}.}
  \item{\dots}{arguments passed to the underlying functions, i.e., to
    \code{\link[partykit]{mob_control}} for \code{raschtree}, and to the
    underlying \code{predict} and \code{plot} methods, respectively.}
  \item{object, x}{an object of class \code{"raschtree"}.}
  \item{newdata}{optional data frame with partitioning variables for
    which predictions should be computed. By default the learning data set
    is used.}
  \item{type}{character specifying the type of predictions or plot. For the
    \code{predict} method, either just the ID of the terminal \code{"node"}
    can be predicted or some property of the model at a given person parameter
    (specified by \code{personpar}).}
  \item{personpar}{numeric person parameter (of length 1) at which the
    predictions are evaluated.}
  \item{terminal_panel, tp_args, tnex, drop_terminal}{arguments passed to
    \code{\link[partykit]{plot.modelparty}}/\code{\link[partykit]{plot.party}}.}
}

\details{
  Rasch trees are an application of model-based recursive partitioning
  (implemented in \code{\link[partykit]{mob}}) to Rasch models
  (implemented in \code{\link[psychotools]{raschmodel}}).
  
  Various methods are provided for \code{"raschtree"} objects, most of them
  inherit their behavior from \code{"modelparty"} objects (e.g., \code{print}, \code{summary},
  etc.). For the Rasch models in the nodes of a tree, \code{coef} extracts all item
  parameters except the first one which is always restricted to be zero. \code{itempar}
  extracts all item parameters (including the first one) and by default restricts their
  sum to be zero (but other restrictions can be used as well). The \code{plot} method
  by default employs the \code{\link{node_profileplot}} panel-generating function and
  the \code{\link{node_regionplot}} panel-generating function is provided as an alternative.

  Rasch tree models are introduced in Strobl et al. (2015), whose analysis
  for the \code{\link{SPISA}} data is replicated in
  \code{vignette("raschtree", package = "psychotree")}. Their illustration
  employing artificial data is replicated below.
}

\references{
  Strobl C, Kopf J, Zeileis A (2015). Rasch Trees:
  A New Method for Detecting Differential Item Functioning in the Rasch
  Model. \emph{Psychometrika}, forthcoming. doi:10.1007/s11336-013-9388-3
}

\value{
  An object of S3 class \code{"raschtree"} inheriting from class \code{"modelparty"}.
}

\seealso{\code{\link[partykit]{mob}}, \code{\link[psychotools]{raschmodel}},
  \code{\link{rstree}}, \code{\link{pctree}}}

\examples{
o <- options(digits = 4)

## artificial data
data("DIFSim", package = "psychotree")

## fit Rasch tree model
rt <- raschtree(resp ~ age + gender + motivation, data = DIFSim)
plot(rt)

## extract item parameters
coef(rt)
itempar(rt)

## inspect parameter stability tests in all splitting nodes
library("strucchange")
sctest(rt, node = 1)
sctest(rt, node = 2)

## highlight items 3 and 14 with DIF
ix <- rep(1, 20)
ix[c(3, 14)] <- 2
plot(rt, ylines = 2.5,  cex = c(0.4, 0.8)[ix], 
  pch = c(19, 19)[ix], col = gray(c(0.5, 0))[ix])

options(digits = o$digits)
}

\keyword{tree}
