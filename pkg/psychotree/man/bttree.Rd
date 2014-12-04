\name{bttree}

\alias{bttree}

\alias{plot.bttree}
\alias{print.bttree}
\alias{predict.bttree}
\alias{itempar.bttree}

\title{Bradley-Terry Tree Models}

\description{
  Recursive partitioning based on Bradley-Terry models.
}

\usage{
bttree(formula, data, na.action,
  type = "loglin", ref = NULL, undecided = NULL, position = NULL, \dots)

\method{predict}{bttree}(object, newdata = NULL,
  type = c("worth", "rank", "best", "node"), \dots)
}

\arguments{
  \item{formula}{A symbolic description of the model to be fit. This
    should be of type \code{y ~ x1 + x2}
    where \code{y} should be an object of class \code{\link[psychotools]{paircomp}}
    and \code{x1} and \code{x2} are used as partitioning variables.}
  \item{data}{an optional data frame containing the variables in the model.}
  \item{na.action}{A function which indicates what should happen when the data
    contain \code{NA}s, defaulting to \code{\link[stats]{na.pass}}.}
  \item{type}{character indicating the type of auxiliary model in \code{bttree} and the
    type of predictions in the \code{predict} method, respectively. For the auxiliary
    model see \code{\link[psychotools]{btmodel}}. For the \code{predict} method, four options are
    available: the fitted \code{"worth"} for each alternative, the corresponding
    \code{"rank"}, the \code{"best"} alternative or the predicted \code{"node"} number.}
  \item{ref, undecided, position}{arguments for the Bradley-Terry
    model passed on to \code{\link[psychotools]{btmodel}}.}
  \item{\dots}{arguments passed to \code{\link[partykit]{mob_control}}.}
  \item{object}{fitted model object of class \code{"bttree"}.}
  \item{newdata}{optionally, a data frame in which to look for variables with
    which to predict. If omitted, the original observations are used.}
}

\details{
  Bradley-Terry tree models are an application of model-based recursive partitioning
  (implemented in \code{\link[partykit]{mob}}) to Bradley-Terry models for
  paired comparison data (implemented in \code{\link[psychotools]{btmodel}}).
  Details about the underlying theory and further explanations of the
  illustrations in the example section can be found in
  Strobl, Wickelmaier, Zeileis (2011).
  
  Various methods are provided for \code{"bttree"} objects, most of them
  inherit their behavior from \code{"mob"} objects (e.g., \code{print}, \code{summary},
  etc.). \code{itempar} behaves analogously to \code{coef} and extracts the
  worth/item parameters from the BT models in the nodes of the tree. The \code{plot}
  method employs the \code{\link{node_btplot}} panel-generating function.
}

\value{
  An object of S3 class \code{"bttree"} inheriting from class \code{"modelparty"}.
}

\references{ 
  Strobl C, Wickelmaier F, Zeileis A (2011). Accounting for
  Individual Differences in Bradley-Terry Models by Means of Recursive Partitioning.
  \emph{Journal of Educational and Behavioral Statistics}, \bold{36}(2), 135--153.
  doi:10.3102/1076998609359791
}

\seealso{\code{\link[partykit]{mob}}, \code{\link[psychotools]{btmodel}}}

\examples{
o <- options(digits = 4)

## Germany's Next Topmodel 2007 data
data("Topmodel2007", package = "psychotree")

## BT tree
tm_tree <- bttree(preference ~ gender + age + q1 + q2 + a3, data = Topmodel2007,
  minsize = 5, ref = "Barbara")
plot(tm_tree, abbreviate = 1, yscale = c(0, 0.5))

## parameter instability tests in root node
library("strucchange")
sctest(tm_tree, node = 1)

## worth parameters in terminal nodes
itempar(tm_tree)

## CEMS university choice data
data("CEMSChoice", package = "psychotree")
summary(CEMSChoice$preference)

## BT tree
cems_tree <- bttree(preference ~ french + spanish + italian + study + work + gender + intdegree,
  data = CEMSChoice, minsize = 5, ref = "London")
plot(cems_tree, abbreviate = 1, yscale = c(0, 0.5))
itempar(cems_tree)

options(digits = o$digits)
}

\keyword{tree}
