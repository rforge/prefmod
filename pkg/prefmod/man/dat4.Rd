\name{dat4}
\alias{dat4}
\docType{data}
\title{Data (paired comparisons): dat4}
\description{
  A fictitious dataset with 100 observations on 6 paired comparisons.
  The responses get the value 1 if the first object in a comparison is preferred and -1 otherwise.                                                                                                     
  For the arrangement of objects and comparisons see \code{\link{llbt.design}}.
}
\usage{data(dat4)}
\format{
  A data frame with 100 observations on 6 comparisons (\code{comp1} to \code{comp6})
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
