\name{bttree}

\alias{bttree}

\alias{plot.bttree}
\alias{print.bttree}
\alias{summary.bttree}
\alias{sctest.bttree}

\alias{deviance.bttree}
\alias{logLik.bttree}
\alias{coef.bttree}
\alias{worth.bttree}
\alias{weights.bttree}

\encoding{latin1}

\title{Bradley-Terry Tree Models}

\description{
  Recursive partitioning based on Bradley-Terry models.
}

\usage{
bttree(formula, data, na.action = na.pass,
  type = "loglin", ref = NULL, undecided = NULL, position = NULL,
  minsplit = 10, \dots)
}

\arguments{
  \item{formula}{A symbolic description of the model to be fit. This
    should be of type \code{y ~ x1 + x2}
    where \code{y} should be an object of class \code{\link[prefmod2]{paircomp}}
    and \code{x1} and \code{x2} are used as partitioning variables.}
  \item{data}{an optional data frame containing the variables in the model.}
  \item{na.action}{A function which indicates what should happen when the data
    contain \code{NA}s, defaulting to \code{\link[stats]{na.pass}}.}
  \item{type, ref, undecided, position}{arguments for the Bradley-Terry
    model passed on to \code{\link[prefmod2]{btreg}}.}
  \item{minsplit, \dots}{arguments passed to \code{\link[party]{mob_control}}.}
}

\details{
  Bradley-Terry tree models are an application of model-based recursive partitioning
  (implemented in \code{\link[party]{mob}}) to Bradley-Terry models for
  paired comparison data (implemented in \code{\link[prefmod2]{btreg}}).
  Internally, a wrapper function \code{\link{btReg}} for \code{btreg} is used.
  For all details about the underlying theory and further explanations
  of the illustrations from the example section can be found in
  Strobl, Wickelmaier, Zeileis (2009).
  
  Various methods are provided for \code{"bttree"} objects, most of them
  inherit their behavior from \code{"mob"} objects (e.g., \code{print}, \code{summary},
  etc.). \code{worth} behaves analogously to \code{coef} and extracts the
  worth parameters from the BT models in the nodes of the tree. The \code{plot}
  method employs the \code{\link{node_btplot}} panel-generating function.
}

\value{
  An object of S3 class \code{"bttree"} which is a list containing only
  a single element of S4 class \code{"mob"} (because this is currently not
  exported).
}

\references{ 
  Carolin Strobl, Florian Wickelmaier, Achim Zeileis (2009). 
  Accounting for Individual Differences in Bradley-Terry Models by Means of Recursive Partitioning.
  Technical Report 54, Department of Statistics, Ludwig-Maximilians-Universität München.
  \url{http://epub.ub.uni-muenchen.de/10588/}
}

\seealso{\code{\link[party]{mob}}, \code{\link[prefmod2]{btreg.fit}}, \code{\link{btReg}}}

\examples{
## package
library("psychotree")

## Germany's Next Topmodel 2007 data
data("Topmodel2007", package = "prefmod2")
## BT tree
tm_tree <- bttree(preference ~ ., data = Topmodel2007, minsplit = 5, ref = "Barbara")
plot(tm_tree, abbreviate = 1, yscale = c(0, 0.5))
## parameter instability tests in root node
sctest(tm_tree, node = 1)
## worth parameters in terminal nodes
worth(tm_tree)

## CEMS university choice data
data("CEMSChoice", package = "prefmod2")
summary(CEMSChoice$preference)
## BT tree
cems_tree <- bttree(preference ~ french + spanish + italian + study + work + gender + intdegree,
  data = CEMSChoice, minsplit = 5, ref = "London")
plot(cems_tree, abbreviate = 1, yscale = c(0, 0.5))
worth(cems_tree)
}

\keyword{tree}
