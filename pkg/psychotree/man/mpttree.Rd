\name{mpttree}

\alias{mpttree}

\alias{plot.mpttree}
\alias{print.mpttree}
\alias{predict.mpttree}

\title{MPT Tree Models}

\description{
  Recursive partitioning based on multinomial processing tree (MPT)models.
}

\usage{
mpttree(formula, data, na.action,
  spec, treeid = NULL,
  optimargs = list(control = list(reltol = .Machine$double.eps^(1/1.2),
                                  maxit = 1000)), \dots)

% \method{predict}{bttree}(object, newdata = NULL,
%   type = c("worth", "rank", "best", "node"), \dots)
}

\arguments{
  \item{formula}{a symbolic description of the model to be fit. This should be
    of type \code{y ~ x1 + x2} where \code{y} should be a matrix of response
    frequencies and \code{x1} and \code{x2} are used as partitioning
    variables.}
  \item{data}{an optional data frame containing the variables in the model.}
  \item{na.action}{a function which indicates what should happen when the data
    contain \code{NA}s, defaulting to \code{\link[stats]{na.pass}}.}
  \item{spec, treeid, optimargs}{arguments for the MPT model passed on to
    \code{\link[psychotools]{mptmodel}}.}
  \item{\dots}{arguments passed to \code{\link[partykit]{mob_control}}.}
% \item{object}{fitted model object of class \code{"bttree"}.}
% \item{newdata}{optionally, a data frame in which to look for variables with
%   which to predict. If omitted, the original observations are used.}
}

\details{
  MPT tree models are an application of model-based recursive partitioning
  (implemented in \code{\link[partykit]{mob}}) to MPT models (implemented in
  \code{\link[psychotools]{mptmodel}}).
  
  Various methods are provided for \code{"mpttree"} objects, most of them
  inherit their behavior from \code{"mob"} objects (e.g., \code{print},
  \code{summary}, etc.).  The \code{plot} method employs the
  \code{\link{node_mptplot}} panel-generating function.
}

\value{
  An object of S3 class \code{"mpttree"} inheriting from class
  \code{"modelparty"}.
}

\references{ 
  Wickelmaier, F., & Zeileis, A. (2013). Investigating the effects of subject
  covariates in multinomial processing tree models by recursive partitioning.
  Presented at the 46th Annual Meeting of the Society of Mathematical
  Psychology, August 4-7, Potsdam, Germany.
}

\seealso{\code{\link[partykit]{mob}}, \code{\link[psychotools]{mptmodel}},
  the \pkg{mpt} package.}

\examples{
o <- options(digits = 4)

## FIX ME
# ## Germany's Next Topmodel 2007 data
# data("Topmodel2007", package = "psychotree")
# 
# ## BT tree
# tm_tree <- bttree(preference ~ ., data = Topmodel2007, minsize = 5, ref = "Barbara")
# plot(tm_tree, abbreviate = 1, yscale = c(0, 0.5))
# 
# ## parameter instability tests in root node
# library("strucchange")
# sctest(tm_tree, node = 1)
# 
# ## worth parameters in terminal nodes
# itempar(tm_tree)
# 
# ## CEMS university choice data
# data("CEMSChoice", package = "psychotree")
# summary(CEMSChoice$preference)
# 
# ## BT tree
# cems_tree <- bttree(preference ~ french + spanish + italian + study + work + gender + intdegree,
#   data = CEMSChoice, minsize = 5, ref = "London")
# plot(cems_tree, abbreviate = 1, yscale = c(0, 0.5))
# itempar(cems_tree)

options(digits = o$digits)
}

\keyword{tree}
