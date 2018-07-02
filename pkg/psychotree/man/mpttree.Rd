\name{mpttree}

\alias{mpttree}

\alias{coef.mpttree}
\alias{plot.mpttree}
\alias{print.mpttree}
\alias{predict.mpttree}

\encoding{latin1}

\title{MPT Trees}

\description{
  Recursive partitioning (also known as trees) based on multinomial processing tree (MPT) models.
}

\usage{
mpttree(formula, data, na.action, cluster, spec, treeid = NULL,
  optimargs = list(control = list(reltol = .Machine$double.eps^(1/1.2),
                                  maxit = 1000)), \dots)
}

\arguments{
  \item{formula}{a symbolic description of the model to be fit. This should be
    of type \code{y ~ x1 + x2} where \code{y} should be a matrix of response
    frequencies and \code{x1} and \code{x2} are used as partitioning
    variables.}
  \item{data}{an optional data frame containing the variables in the model.}
  \item{na.action}{a function which indicates what should happen when the data
    contain \code{NA}s, defaulting to \code{\link[stats]{na.pass}}.}
  \item{cluster}{optional vector (typically numeric or factor) with a
    cluster ID to be employed for clustered covariances in the parameter
    stability tests.}
  \item{spec, treeid, optimargs}{arguments for the MPT model passed on to
    \code{\link[psychotools]{mptmodel}}.}
  \item{\dots}{arguments passed to \code{\link[partykit]{mob_control}}.}
}

\details{
  MPT trees (Wickelmaier & Zeileis, 2018) are an application of
  model-based recursive partitioning (implemented in
  \code{\link[partykit]{mob}}) to MPT models (implemented in
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
Wickelmaier F, Zeileis A (2018). Using Recursive Partitioning to Account
  for Parameter Heterogeneity in Multinomial Processing Tree Models.
  \emph{Behavior Research Methods}, \bold{50}(3), 1217--1233.
  \doi{10.3758/s13428-017-0937-z}
}

\seealso{\code{\link[partykit]{mob}}, \code{\link[psychotools]{mptmodel}}.}

\examples{
o <- options(digits = 4)

## Source Monitoring data
data("SourceMonitoring", package = "psychotools")

## MPT tree
sm_tree <- mpttree(y ~ sources + gender + age, data = SourceMonitoring,
  spec = mptspec("SourceMon", .restr = list(d1 = d, d2 = d)))
plot(sm_tree, index = c("D1", "D2", "d", "b", "g"))

## extract parameter estimates
coef(sm_tree)

## parameter instability tests in root node
library("strucchange")
sctest(sm_tree, node = 1)

## storage and retrieval deficits in psychiatric patients
data("MemoryDeficits", package = "psychotools")
MemoryDeficits$trial <- ordered(MemoryDeficits$trial)

## MPT tree
sr_tree <- mpttree(cbind(E1, E2, E3, E4) ~ trial + group,
  data = MemoryDeficits, cluster = ID, spec = mptspec("SR2"), alpha = 0.1)

## extract parameter estimates
coef(sr_tree)

options(digits = o$digits)
}

\keyword{tree}
