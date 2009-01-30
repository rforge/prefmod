\name{dat4}
\alias{dat4}
\docType{data}
\title{Data (paired comparisons): dat4}
\description{
  A fictitious data set with 100 observations on 6 paired comparisons.
  The responses are 1 assigned if the first object in a comparison is the preferred one and -1 if it is not preferred.
  For the ordering of objects and comparisons see \code{\link{llbt.design}}.
}
\usage{data(dat4)}
\format{
  A data frame with 100 observations on 6 comparisons:
  \describe{
    \item{\code{comp1}}{responses to first comparison}
    \item{\code{comp2}}{responses to second comparison}
    \item{\code{comp3}}{responses to third comparison}
    \item{\code{comp4}}{responses to fourth comparison}
    \item{\code{comp5}}{responses to fifth comparison}
    \item{\code{comp6}}{responses to sixth comparison}
  }
}
\examples{
data(dat4)
str(dat4)
# to get a general idea we use the histogram plot
par(mfrow=c(2,3))
for (i in 1:6) barplot(table(dat4[,i]))
par(mfrow=c(1,1))
}
\keyword{datasets}
