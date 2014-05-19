\encoding{UTF-8}
\name{salad}
\alias{salad}
\docType{data}
\title{Data (ranks): Salad Dressings (Critchlow and Fligner)}
\description{
  The dataset contains the rankings of four salad dressings concerning tartness by 32 judges, 
  with values ranging from 1 (most tart) 
  to 4 (least tart).
}
\usage{salad}
\format{
  A data frame with 32 observations on 4 variables (\code{A}, \code{B}, \code{C}, \code{D}) each
 representing a different salad dressing.
}
\references{
  Douglas E. Crichtlow and Michael A. Fligner (1991), Paired Comparison, Triple Comparison, and Ranking Expermients As Generalized Linear Models, and Their Implementation on GLIM. Psychometrika Vol. 56, No. 3
}
\examples{
# Example for object covariates
# fit object covariates:
# salads A - D have varying concentrations of acetic and gluconic acid.
# The four pairs of concentrations are
# A = (.5, 0), B = (.5, 10.0), C = (1.0, 0), and D = (0, 10.0),

conc <- matrix(c(.5,0, .5,10, 1,0, 0,10), ncol = 2, byrow = TRUE)
sal  <- patt.design(salad, nitems = 4, resptype = "ranking")
X    <- as.matrix(sal[, 2:5]) %*% conc

glm(y ~ X, data = sal, family = poisson)
}
\keyword{datasets}
