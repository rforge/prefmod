\name{rstree}
\alias{rstree}

\alias{plot.rstree}
\alias{print.rstree}
\alias{predict.rstree}
\alias{itempar.rstree}
\alias{threshpar.rstree}

\title{Rating Scale Tree Models}

\description{
  Recursive partitioning based on rating scale models.
}

\usage{
rstree(formula, data, na.action, reltol = 1e-10,
  deriv = c("sum", "diff"), maxit = 100L, \dots)

\method{predict}{rstree}(object, newdata = NULL,
  type = c("probability", "cumprobability", "mode", "median", "mean",
    "category-information", "item-information", "test-information", "node"),
  personpar = 0, \dots)

\method{plot}{rstree}(x, type = c("regions", "profile"), terminal_panel = NULL,
  tp_args = list(...), tnex = 2L, drop_terminal = TRUE, ...)
}

\arguments{
  \item{formula}{A symbolic description of the model to be fit. This
    should be of type \code{y ~ x1 + x2} where \code{y} should be a
    matrix with items in the columns and observations in the rows and
    \code{x1} and \code{x2} are used as partitioning
    variables. Additionally each item (column) should have the same
    maximum value (see \code{\link{pctree}} for a way to handle
    variable maximum values).}
  \item{data}{a data frame containing the variables in the model.}
  \item{na.action}{a function which indicates what should happen when the data
    contain missing values (\code{NA}s).}
  \item{deriv}{character. If "sum" (the default), the first derivatives
    of the elementary symmetric functions are calculated with the sum
    algorithm. Otherwise ("diff") the difference algorithm (faster but
    numerically unstable) is used.}
  \item{reltol, maxit}{arguments passed via \code{\link[psychotools]{rsmodel}}
    to \code{\link[stats]{optim}}.}
  \item{\dots}{arguments passed to the underlying functions, i.e., to
    \code{\link[partykit]{mob_control}} for \code{rstree}, and to the
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
  Rating scale tree models are an application of model-based recursive partitioning
  (implemented in \code{\link[partykit]{mob}}) to rating scale models
  (implemented in \code{\link[psychotools]{rsmodel}}).
  
  Various methods are provided for \code{"rstree"} objects, most of them
  inherit their behavior from \code{"mob"} objects (e.g., \code{print},
  \code{summary}, etc.). For the rating scale models in the nodes of a
  tree, \code{coef} extracts all item parameters. The \code{plot} method
  employs the \code{\link{node_regionplot}} panel-generating function by
  default.

  Various methods are provided for \code{"rstree"} objects, most of them
  inherit their behavior from \code{"modelparty"} objects (e.g., \code{print}, \code{summary},
  etc.). For the RSMs in the nodes of a tree, \code{coef} extracts all item and threshold
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
  An object of S3 class \code{"rstree"} inheriting from class \code{"modelparty"}.
}

\seealso{\code{\link[partykit]{mob}}, \code{\link[psychotools]{rsmodel}},
  \code{\link{pctree}}, \code{\link{raschtree}}}

\examples{
## verbal aggression data from package psychotools
data("VerbalAggression", package = "psychotools")

## responses to the first other-to-blame situation (bus)
VerbalAggression$s1 <- VerbalAggression$resp[, 1:6]

## exclude subjects who only scored in the highest or the lowest categories
VerbalAggression <- subset(VerbalAggression, rowSums(s1) > 0 & rowSums(s1) < 12)

## fit rating scale tree model for the first other-to-blame situation
rst <- rstree(s1 ~ anger + gender, data = VerbalAggression)

## print tree (with and without parameters)
print(rst)
print(rst, FUN = function(x) " *")

## show summary for terminal panel nodes
summary(rst)

## visualization
plot(rst, type = "regions")
plot(rst, type = "profile")

## extract item and threshold parameters
coef(rst)
itempar(rst)
threshpar(rst)

## inspect parameter stability tests in all splitting nodes
if(require("strucchange")) {
sctest(rst, node = 1)
sctest(rst, node = 2)
}
}

\keyword{tree}
