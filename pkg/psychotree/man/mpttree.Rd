\name{mpttree}
\alias{mpttree}

\alias{plot.mpttree}
\alias{print.mpttree}
\alias{summary.mpttree}
\alias{sctest.mpttree}

\alias{deviance.mpttree}
\alias{logLik.mpttree}
\alias{coef.mpttree}
\alias{weights.mpttree}

\title{MPT Tree Models}

\description{
  Recursive partitioning based on multinomial processing tree (MPT) models.
}

\usage{
mpttree(formula, data, treeid = NULL, na.action = na.pass, minsplit = 5,
  \dots)
}

\arguments{
  \item{formula}{A symbolic description of the model to be fit. This
    should be of type \code{y ~ list(eq1, eq2, ...) | x1 + x2}
    where \code{y} should be a matrix of response frequencies, \code{eqX} is
    a symbolic representaion of the probability of a response in category
    \code{X} and \code{x1} and \code{x2} are used as partitioning variables.}
  \item{data}{a data frame containing the variables in the model.}
  \item{treeid}{currently not used.}
  \item{na.action}{currently not used.}
  \item{minsplit, \dots}{arguments passed to \code{\link[party]{mob_control}}.}
}

\details{
  MPT tree models are an application of model-based recursive partitioning
  (implemented in \code{\link[party]{mob}}) to multinomial processing tree
  (MPT) models (implemented in \code{\link{mptModel}}).
  
  Various methods are provided for \code{"mpttree"} objects, most of them
  inherit their behavior from \code{"mob"} objects (e.g., \code{print},
  \code{summary}, etc.). For the MPT models in the nodes of a tree,
  \code{coef} extracts all model parameters. The \code{plot} method employs
  the \code{\link{node_mptplot}} panel-generating function.

% Rasch tree models are introduced in Strobl et al. (2010), whose analysis
% for the \code{\link{SPISA}} data is replicated in
% \code{vignette("raschtree", package = "psychotree")}. Their illustration
% employing artificial data is replicated below.
}

% \references{
%   Strobl, C., Kopf, J., and Zeileis, A. (2010). A New Method for Detecting
%   Differential Item Functioning in the Rasch Model. Technical Report 92.
%   \url{http://epub.ub.uni-muenchen.de/11915/}
% }

\value{
  An object of S3 class \code{"mpttree"} which is a list containing only
  a single element of S4 class \code{"mob"} (because this is currently not
  exported from the party package).
}

\seealso{\code{\link[party]{mob}}, \code{\link{mptModel}}, \code{\link{bttree}}}

\examples{
## TODO
# ## artificial data
# data("DIFSim", package = "psychotree")
# 
# ## fit Rasch tree model
# rt <- raschtree(resp ~ age + gender + motivation, data = DIFSim)
# plot(rt)
# 
# ## extract item parameters
# coef(rt)
# worth(rt)
# 
# ## inspect parameter stability tests in all splitting nodes
# sctest(rt, node = 1)
# sctest(rt, node = 2)
# 
# ## highlight items 3 and 14 with DIF
# ix <- rep(1, 20)
# ix[c(3, 14)] <- 2
# plot(rt, ylines = 2.5,  cex = c(0.4, 0.8)[ix], 
#   pch = c(19, 19)[ix], col = gray(c(0.5, 0))[ix])
}

\keyword{tree}
