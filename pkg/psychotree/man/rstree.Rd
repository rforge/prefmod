\name{rstree}

\alias{rstree}

\alias{plot.rstree}
\alias{print.rstree}
\alias{summary.rstree}
\alias{sctest.rstree}
\alias{deviance.rstree}
\alias{logLik.rstree}
\alias{coef.rstree}
\alias{weights.rstree}


\title{Rating Scale Tree Models}


\description{

  Recursive partitioning based on rating scale models.

}


\usage{

rstree(formula, data, minsplit = 20, gradtol = 1e-6,
  deriv = c("sum", "diff"), hessian = TRUE, iterlim = 100L, \dots)

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

  \item{deriv}{character. If "sum" (the default), the first derivatives
    of the elementary symmetric functions are calculated with the sum
    algorithm. Otherwise ("diff") the difference algorithm (faster but
    numerically unstable) is used.}
  
  \item{hessian}{logical. Should the Hessian of the final model be computed?
    If set to \code{FALSE}, the \code{vcov} method can only return \code{NA}s
    and consequently no standard errors or tests are available in the \code{summary}.}

  \item{gradtol, iterlim}{arguments passed via \code{\link{RSModel}}
    to \code{\link[stats]{nlm}}.}

  \item{minsplit, \dots}{arguments passed to
  \code{\link[party]{mob_control}}.}

}


\details{

  Rating scale tree models are an application of model-based recursive partitioning
  (implemented in \code{\link[party]{mob}}) to rating scale models
  (implemented in \code{\link{RSModel}}).
  
  Various methods are provided for \code{"rstree"} objects, most of them
  inherit their behavior from \code{"mob"} objects (e.g., \code{print},
  \code{summary}, etc.). For the rating scale models in the nodes of a
  tree, \code{coef} extracts all item parameters. The \code{plot} method
  employs the \code{\link{node_effects}} panel-generating function by
  default.

  %% Alternatively, this can be also set to
  %% \code{\link{node_effects}} or \code{\link{node_ccc}} via argument
  %% \code{terminal_panel} of the plot method. (See
  %% \code{\link[party]{plot.BinaryTree}} for details).

  % Rasch tree models are introduced in Strobl et al. (2010), whose analysis
  % for the \code{\link{SPISA}} data is replicated in
  % \code{vignette("rstree", package = "psychotree")}. Their illustration
  % employing artificial data is replicated below.
}


\references{
  Abou El-Komboz, B., Zeileis, A., and Strobl, C. (2014). Detecting
  Differential Item and Step Functioning with Rating Scale and Partial
  Credit Trees. Technical Report 152. Department of Statistics, Ludwig-
  Maximilians-Universitaet Muenchen. \url{http://epub.ub.uni-muenchen.de/17984/}
}


\value{

  An object of S3 class \code{c("rstree", "pctree")} which is a list containing only
  a single element of S4 class \code{"mob"} (because this is currently not
  exported from the party package).

}


\seealso{\code{\link[party]{mob}}, \code{\link{RSModel}},
  \code{\link{pctree}}, \code{\link{raschtree}}}


\examples{

## verbal aggression data from package psychotools
data("VerbalAggression", package = "psychotools")

## exclude subjects who only scored in the highest or the lowest categories
VerbalAggression <- subset(VerbalAggression,
                             rowSums(resp[, 1:6]) > 0 & rowSums(resp[, 1:6]) < 12)

## fit rating scale tree model for the first other-to-blame situation
rst <- rstree(resp[, 1:6, drop = FALSE] ~ anger + gender, data =
VerbalAggression)

## just print mob object
print(rst)

## show summary for terminal panel nodes
summary(rst)

## plot with default (and currently only) terminal panel function (node_effects)
plot(rst)

## extract item and category parameters
coef(rst)

## inspect parameter stability tests in all splitting nodes
sctest(rst, node = 1)
sctest(rst, node = 2)

}


\keyword{tree}
