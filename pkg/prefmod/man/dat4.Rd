\encoding{UTF-8}
\name{dat4}
\alias{dat4}
\docType{data}
\title{Data (paired comparisons): dat4}
\description{%
  A fictitious dataset with 100 observations on 6 paired comparisons.
  The responses get the value \eqn{1} if the first object in a comparison is preferred and \eqn{-1} otherwise.

  For the arrangement of objects and comparisons see \code{\link{llbt.design}}.
}
\usage{dat4}
\format{%
  A data frame with 100 observations on 6 comparisons (\code{comp1} to \code{comp6})
}
\examples{
str(dat4)

# to get a general idea we use the histogram plot
old_par <- par(mfrow = c(2, 3))
for(i in 1:6){ barplot(table(dat4[, i])) }
par(old_par)
}
\keyword{datasets}
