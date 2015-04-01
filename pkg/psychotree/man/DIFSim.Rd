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
  and Abou El-Komboz et al. (2014). See the manual pages for
  \code{\link{raschtree}} and \code{\link{pctree}} for fitting the
  tree models..
}

\references{
  Abou El-Komboz B, Zeileis A, Strobl C (2014). Detecting
  Differential Item and Step Functioning with Rating Scale and Partial
  Credit Trees. Technical Report 152. Department of Statistics,
  Ludwig-Maximilians-Universitaet Muenchen.
  \url{http://epub.ub.uni-muenchen.de/17984/}

  Strobl C, Kopf J, Zeileis A (2015). Rasch Trees:
  A New Method for Detecting Differential Item Functioning in the Rasch
  Model. \emph{Psychometrika}, forthcoming. doi:10.1007/s11336-013-9388-3
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
