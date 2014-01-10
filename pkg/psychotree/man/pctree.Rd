\name{pctree}

\alias{pctree}

\alias{plot.pctree}
\alias{print.pctree}
\alias{summary.pctree}
\alias{sctest.pctree}
\alias{deviance.pctree}
\alias{logLik.pctree}
\alias{coef.pctree}
\alias{weights.pctree}


\title{Partial Credit Tree Models}


\description{

  Recursive partitioning based on partial credit models.

}


\usage{

pctree(formula, data, minsplit = 30, diff = FALSE,
  nullcats = c("keep", "downcode", "ignore"), start = NULL,
  gradtol = 1e-6, iterlim = 100L, hessian = TRUE, \dots)

}


\arguments{

  \item{formula}{A symbolic description of the model to be fit. This
    should be of type \code{y ~ x1 + x2} where \code{y} should be a
    matrix with items in the columns and observations in the rows and
    \code{x1} and \code{x2} are used as partitioning
    variables.}

  \item{data}{a data frame containing the variables in the model.}

  \item{diff}{logical. If \code{FALSE}, the first derivatives of the
    elementary symmetric functions are calculated with the sum
    algorithm. Otherwise (\code{TRUE}) the difference algorithm is used.}

  \item{nullcats}{character. How null categories should be treated. See
  \code{\link[psychotools]{PCModel.fit}} for details.}
  
  \item{start, gradtol, iterlim}{arguments passed via \code{\link{PCModel}}
    to \code{\link[stats]{nlm}}.}
  
  \item{hessian}{logical. Should the Hessian of the final model be computed?
    If set to \code{FALSE}, the \code{vcov} method can only return \code{NA}s
    and consequently no standard errors or tests are available in the \code{summary}.}

  \item{minsplit, \dots}{arguments passed to \code{\link[party]{mob_control}}.}

}


\details{

  Partial credit tree models are an application of model-based recursive partitioning
  (implemented in \code{\link[party]{mob}}) to partial credit models
  (implemented in \code{\link{PCModel}}).
  
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


% \references{
%   Strobl, C., Kopf, J., and Zeileis, A. (2010). A New Method for Detecting
%   Differential Item Functioning in the Rasch Model. Technical Report 92.
%   Department of Statistics, Ludwig-Maximilians-Universitaet Muenchen.
%   \url{http://epub.ub.uni-muenchen.de/11915/}
% }


\value{

  An object of S3 class \code{"pctree"} which is a list containing only
  a single element of S4 class \code{"mob"} (because this is currently not
  exported from the party package).

}


\seealso{\code{\link[party]{mob}}, \code{\link{PCModel}},
  \code{\link{rstree}}, \code{\link{raschtree}}}


\examples{

## verbal aggression data from package psychotools
data("VerbalAggression", package = "psychotools")

## exclude subjects who only scored in the highest or the lowest categories
VerbalAggression <- subset(VerbalAggression,
                            rowSums(resp[, 7:12]) > 0 & rowSums(resp[, 7:12]) < 12)

## fit partial credit tree model for the second other-to-blame situation
pct <- pctree(resp[, 7:12, drop = FALSE] ~ anger + gender, data =
VerbalAggression)

## just print mob object
print(pct)

## show summary for terminal panel nodes
summary(pct)

## plot with default (and currently only) terminal panel function (node_effects)
plot(pct)

## extract item and category parameters
coef(pct)

## inspect parameter stability tests in the splitting node
sctest(pct, node = 1)

}


\keyword{tree}
