\name{prefmod-package}
\alias{prefmod-package}
\alias{prefmod}
\docType{package}
\title{
Utilities to fit paired comparison models for preferences
}
\description{
Generates design matrix for analysing real paired comparisons and derived
paired comparison data (Likert type items / ratings or rankings) using a
loglinear approach. Fits loglinear Bradley-Terry model (LLBT) exploting an
eliminate feature. Computes pattern models for paired comparisons,
rankings, and ratings. Some treatment of missing values (MCAR and MNAR).
Fits pattern mixture models using a nonparametric ML approach.
}
\details{
\tabular{ll}{
Package: \tab prefmod\cr
Type: \tab Package\cr
Version: \tab 0.8-24\cr
Date: \tab 2011-03-29\cr
Depends: \tab stats\cr
License: \tab GPL\cr
}
}
\author{
Reinhold Hatzinger

Maintainer: Reinhold Hatzinger <reinhold.hatzinger@wu.ac.at>
}
\keyword{ package }
\examples{
## mini example with three Likert items
## and two subject covariates

data(xmpl)   # example data in package
dsgnmat <- patt.design(xmpl, nitems = 3, resptype="rating",
    blnIntcovs = TRUE, cov.sel="ALL")
print(head(dsgnmat))

## fit of Critchlov & Fligner (1991) Salad Dressings Data
data(salad)
pattR.fit(salad, nitems=4)

# alternatively use glm() with patt.design()
sal<-patt.design(salad,nitems=4,resptype="ranking")
glm(y~A+B+C+D,family=poisson,data=sal)
}
