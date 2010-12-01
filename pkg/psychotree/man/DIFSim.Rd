\name{DIFSim}
\alias{DIFSim}
\encoding{latin1}

\title{Artificial Data with Differential Item Functioning}

\description{
  Artificial data simulated from a Rasch model where the items
  exhibit differential item functioning (DIF).
}

\usage{data("DIFSim")}

\format{
A data frame containing 200 observations on 4 variables.
\describe{
  \item{resp}{matrix with \code{0}/\code{1} results for 20 items.}
  \item{age}{age in years.}
  \item{gender}{factor indicating gender.}
  \item{motivation}{ordered factor indicating motivation level.}
}
}

\details{  
  The data are employed for illustration in Strobl et al. (2010),
  whose results are replicated on the manual page for
  \code{\link{raschtree}}.
}

\references{
  Strobl, C., Kopf, J., and Zeileis, A. (2010). A New Method for Detecting
  Differential Item Functioning in the Rasch Model. Technical Report 92.
  Department of Statistics, Ludwig-Maximilians-Universität München.
  \url{http://epub.ub.uni-muenchen.de/11915/}
}

\seealso{\code{\link{raschtree}}}

\examples{
## data
data("DIFSim", package = "psychotree")

## summary of covariates
summary(DIFSim[,-1])

## histogram of raw scores
hist(rowSums(DIFSim$resp), breaks = 0:19 + 0.5)
}

\keyword{datasets}
