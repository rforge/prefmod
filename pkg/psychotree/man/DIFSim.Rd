\name{DIFSim}
\alias{DIFSim}
\alias{DIFSimPC}

\title{Artificial Data with Differential Item Functioning}

\description{
  Artificial data simulated from a Rasch model and a partial credit model,
  respectively, where the items exhibit differential item functioning (DIF).
}

\usage{
data(DIFSim)
data(DIFSimPC)
}

\format{
Two data frames containing 200 and 500 observations, respectively, on 4 variables.
\describe{
  \item{resp}{an \code{\link[psychotools]{itemresp}} matrix with binary or polytomous
    results for 20 or 8 items, respectively.}
  \item{age}{age in years.}
  \item{gender}{factor indicating gender.}
  \item{motivation}{ordered factor indicating motivation level.}
}
}

\details{  
  The data are employed for illustrations in Strobl et al. (2015)
  and Komboz et al. (2018). See the manual pages for
  \code{\link{raschtree}} and \code{\link{pctree}} for fitting the
  tree models..
}

\references{
Komboz B, Zeileis A, Strobl C (2018).
  Tree-Based Global Model Tests for Polytomous Rasch Models.
  \emph{Educational and Psychological Measurement}, \bold{78}(1), 128--166.
  \doi{10.1177/0013164416664394}

Strobl C, Kopf J, Zeileis A (2015). Rasch Trees: A New Method for Detecting
  Differential Item Functioning in the Rasch Model.
  \emph{Psychometrika}, \bold{80}(2), 289--316.
  \doi{10.1007/s11336-013-9388-3}
}

\seealso{\code{\link{raschtree}}, \code{\link{pctree}}}

\examples{
## data
data("DIFSim", package = "psychotree")
data("DIFSimPC", package = "psychotree")

## summary of covariates
summary(DIFSim[, -1])
summary(DIFSimPC[, -1])

## empirical frequencies of responses
plot(DIFSim$resp)
plot(DIFSimPC$resp)

## histogram of raw scores
hist(rowSums(DIFSim$resp), breaks = 0:20 - 0.5)
hist(rowSums(DIFSimPC$resp), breaks = 0:17 - 0.5)
}

\keyword{datasets}
