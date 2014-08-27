\name{pctree}
\alias{pctree}

\alias{plot.pctree}
\alias{print.pctree}
\alias{predict.pctree}
\alias{itempar.pctree}
\alias{threshpar.pctree}

\title{Partial Credit Tree Models}

\description{
  Recursive partitioning based on partial credit models.
}

\usage{
pctree(formula, data, na.action, nullcats = c("keep", "downcode", "ignore"),
  reltol = 1e-10,  deriv = c("sum", "diff"), maxit = 100L, \dots)

\method{predict}{pctree}(object, newdata = NULL,
  type = c("probability", "cumprobability", "mode", "median", "mean",
    "category-information", "item-information", "test-information", "node"),
  personpar = 0, \dots)

\method{plot}{pctree}(x, type = c("regions", "profile"), terminal_panel = NULL,
  tp_args = list(...), tnex = 2L, drop_terminal = TRUE, ...)
}

\arguments{
  \item{formula}{A symbolic description of the model to be fit. This
    should be of type \code{y ~ x1 + x2} where \code{y} should be a
    matrix with items in the columns and observations in the rows and
    \code{x1} and \code{x2} are used as partitioning
    variables.}
  \item{data}{a data frame containing the variables in the model.}
  \item{na.action}{a function which indicates what should happen when the data
    contain missing values (\code{NA}s).}
  \item{nullcats}{character. How null categories should be treated. See
  \code{\link[psychotools]{pcmodel}} for details.}
  \item{deriv}{character. If "sum" (the default), the first derivatives
    of the elementary symmetric functions are calculated with the sum
    algorithm. Otherwise ("diff") the difference algorithm (faster but
    numerically unstable) is used.}
  \item{reltol, maxit}{arguments passed via \code{\link[psychotools]{pcmodel}}
    to \code{\link[stats]{optim}}.}
  \item{\dots}{arguments passed to the underlying functions, i.e., to
    \code{\link[partykit]{mob_control}} for \code{pctree}, and to the
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
  Partial credit tree models are an application of model-based recursive partitioning
  (implemented in \code{\link[partykit]{mob}}) to partial credit models
  (implemented in \code{\link[psychotools]{pcmodel}}).
  
  Various methods are provided for \code{"pctree"} objects, most of them
  inherit their behavior from \code{"modelparty"} objects (e.g., \code{print}, \code{summary},
  etc.). For the PCMs in the nodes of a tree, \code{coef} extracts all item and threshold
  parameters except those restricted to be zero. \code{itempar} and \code{threshpar}
  extract all item and threshold parameters (including the restricted ones).
  The \code{plot} method by default employs the \code{\link{node_regionplot}}
  panel-generating function and the \code{\link{node_profileplot}} panel-generating
  function is provided as an alternative.
}

\references{
  Abou El-Komboz B, Zeileis A, Strobl C (2014). Detecting
  Differential Item and Step Functioning with Rating Scale and Partial
  Credit Trees. Technical Report 152. Department of Statistics,
  Ludwig-Maximilians-Universitaet Muenchen.
  \url{http://epub.ub.uni-muenchen.de/17984/}
}

\value{
  An object of S3 class \code{"pctree"} inheriting from class \code{"modelparty"}.
}

\seealso{\code{\link[partykit]{mob}}, \code{\link[psychotools]{pcmodel}},
  \code{\link{rstree}}, \code{\link{raschtree}}}

\examples{
## verbal aggression data from package psychotools
data("VerbalAggression", package = "psychotools")

## use response to the second other-to-blame situation (train)
VerbalAggression$s2 <- VerbalAggression$resp[, 7:12]

## exclude subjects who only scored in the highest or the lowest categories
VerbalAggression <- subset(VerbalAggression, rowSums(s2) > 0 & rowSums(s2) < 12)

## fit partial credit tree model
pct <- pctree(s2 ~ anger + gender, data = VerbalAggression)

## print tree (with and without parameters)
print(pct)
print(pct, FUN = function(x) " *")

## show summary for terminal panel nodes
summary(pct)

## visualization
plot(pct, type = "regions")
plot(pct, type = "profile")

## extract item and threshold parameters
coef(pct)
itempar(pct)
threshpar(pct)

## inspect parameter stability tests in the splitting node
if(require("strucchange")) {
sctest(pct, node = 1)
}
}

\keyword{tree}
