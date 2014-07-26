\name{pctree}
\alias{pctree}

\alias{plot.pctree}
\alias{print.pctree}

\title{Partial Credit Tree Models}

\description{
  Recursive partitioning based on partial credit models.
}

\usage{
pctree(formula, data, na.action = na.pass, nullcats = c("keep", "downcode", "ignore"),
  reltol = 1e-10,  deriv = c("sum", "diff"), maxit = 100L, \dots)
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
  
  \item{\dots}{arguments passed to \code{\link[partykit]{mob_control}}.}
}


\details{
  Partial credit tree models are an application of model-based recursive partitioning
  (implemented in \code{\link[partykit]{mob}}) to partial credit models
  (implemented in \code{\link[psychotools]{pcmodel}}).
  
  Various methods are provided for \code{"pctree"} objects, most of them
  inherit their behavior from \code{"mob"} objects (e.g., \code{print},
  \code{summary}, etc.). For the partial credit models in the nodes of a
  tree, \code{coef} extracts all threshold parameters. The \code{plot}
  method employs the \code{\link{node_effects}} panel-generating
  function by default.

  %% Alternatively, this can be also set to
  %% \code{\link{node_ccc}} via argument \code{terminal_panel} of the plot
  %% method. (See \code{\link[party]{plot.BinaryTree}} for details).

  % Rasch tree models are introduced in Strobl et al. (2010), whose analysis
  % for the \code{\link{SPISA}} data is replicated in
  % \code{vignette("pctree", package = "psychotree")}. Their illustration
  % employing artificial data is replicated below.
}

\references{
  Abou El-Komboz, B., Zeileis, A., and Strobl, C. (2014). Detecting
  Differential Item and Step Functioning with Rating Scale and Partial
  Credit Trees. Technical Report 152. Department of Statistics, Ludwig-
  Maximilians-Universitaet Muenchen. \url{http://epub.ub.uni-muenchen.de/17984/}
}

\value{
  An object of S3 class \code{"pctree"} which is a list containing only
  a single element of S4 class \code{"mob"} (because this is currently not
  exported from the party package).
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

## just print mob object
print(pct)

## show summary for terminal panel nodes
summary(pct)

## plot with default (and currently only) terminal panel function (node_effects)
plot(pct)

## extract item and category parameters
coef(pct)

## inspect parameter stability tests in the splitting node
library("strucchange")
sctest(pct, node = 1)
}

\keyword{tree}
