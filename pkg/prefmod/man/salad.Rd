\name{salad}
\alias{salad}
\docType{data}
\title{Data (ranks): Salad Dressings (Critchlow and Fligner)}
\description{
  The data frame shows the ranking of four salad dressing concerning tartness by 32 judges, whereas (1) is assigned as most tart and (4) is assigned as least tart.
}
\usage{data(salad)}
\format{
\describe{
    \item{\code{A}}{salad dressing preperation A}
    \item{\code{B}}{salad dressing preperation B}
    \item{\code{C}}{salad dressing preperation C}
    \item{\code{D}}{salad dressing preperation D}
  }
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
data(salad)
conc<-matrix(c(.5,0, .5,10,  1,0, 0,10), nc=2, byrow=TRUE)
sal<-patt.design(salad,nitems=4,resptype="ranking")
X<-as.matrix(sal[,2:5])%*%conc
glm(y~X,family=poisson,data=sal)
}
\keyword{datasets}
