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
mpttree(formula, data, na.action = na.pass, minsplit = 5, maxit = 1000,
        \dots)
}

\arguments{
  \item{formula}{A symbolic description of the model to be fit. This
    should be of type \code{y ~ list(eq1, eq2, ...) | x1 + x2}
    where \code{y} should be a matrix of response frequencies, \code{eqX} is
    a symbolic representaion of the probability of a response in category
    \code{X} and \code{x1} and \code{x2} are used as partitioning variables.}
  \item{data}{a data frame containing the variables in the model.}
  \item{na.action}{currently not used.}
  \item{minsplit, \dots}{arguments passed to \code{\link[party]{mob_control}}.}
  \item{maxit}{the maximum number of iterations, passed to
    \code{\link[psychotools]{mptModel.fit}}.}
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

\references{
  Wickelmaier, F., and Zeileis, A. (2013). A first implementation of
  recursive partitioning for multinomial processing tree models. Presented
  at the Psychoco 2013 International Workshop on Psychometric Computing,
  February 14-15, Zurich, Switzerland. 
}

\value{
  An object of S3 class \code{"mpttree"} which is a list containing only
  a single element of S4 class \code{"mob"} (because this is currently not
  exported from the party package).
}

\seealso{\code{\link[party]{mob}}, \code{\link{mptModel}}, \code{\link{bttree}}}

\examples{
## source monitoring performance for different sources
data("SourceMonitoring", package = "psychotree")

## MPT tree
sm_tree <- mpttree(cbind(aa, ab, an, ba, bb, bn, na, nb, nn) ~
  list(
    D1*d + D1*(1 - d)*g + (1 - D1)*b*g,
    D1*(1 - d)*(1 - g) + (1 - D1)*b*(1 - g),
    (1 - D1)*(1 - b),
    D2*(1 - d)*g + (1 - D2)*b*g,
    D2*d + D2*(1 - d)*(1 - g) + (1 - D2)*b*(1 - g),
    (1 - D2)*(1 - b),
    b*g,
    b*(1 - g),
    1 - b
) | sources + age + gender, SourceMonitoring)
plot(sm_tree)

## extract parameter estimates
coef(sm_tree)

## inspect parameter stability tests in all splitting nodes
sctest(sm_tree, node = 1)
sctest(sm_tree, node = 2)

## storage and retrieval deficits in psychiatric patients
data("MemoryDeficits", package = "psychotree")

## MPT tree
sr_tree <- mpttree(cbind(E1, E2, E3, E4) ~
  list(
    c*r,
    (1 - c)*u^2,
    2*(1 - c)*u*(1 - u),
    c*(1 - r) + (1 - c)*(1 - u)^2
) | trial + group, MemoryDeficits, maxit = 3000)

## extract parameter estimates
coef(sr_tree)
}

\keyword{tree}
