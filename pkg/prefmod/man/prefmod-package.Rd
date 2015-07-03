\encoding{UTF-8}
\name{prefmod-package}
\alias{prefmod-package}
\alias{prefmod}
\docType{package}
\title{\pkg{prefmod}: Utilities to Fit Paired Comparison Models for Preferences}
\description{
Generates design matrix for analysing real paired comparisons and derived
paired comparison data (Likert type items / ratings or rankings) using a
loglinear approach. Fits loglinear Bradley-Terry model (LLBT) exploiting an
eliminate feature. Computes pattern models for paired comparisons,
rankings, and ratings. Some treatment of missing values (MCAR and MNAR).
Fits pattern mixture models using a non-parametric ML approach.
}
\details{\tabular{ll}{
Package:\tab{}prefmod\cr%
Type:   \tab{}Package\cr%
Version:\tab{}\Sexpr[stage=build]{packageDescription("prefmod")$Version}\cr%
Date:   \tab{}\Sexpr[stage=build]{packageDescription("prefmod")$Date}\cr%
License:\tab{}\Sexpr[stage=build]{packageDescription("prefmod")$License}\cr%
}}
\author{Reinhold Hatzinger, Marco J. Maier

Maintainer: Marco J. Maier (\email{marco.maier@wu.ac.at})%
}
\references{Hatzinger, R., & Dittrich, R. (2012). \pkg{prefmod}: An R Package for Modeling Preferences Based on Paired Comparisons, Rankings, or Ratings.
\emph{Journal of Statistical Software, 48}(10), 1--31. \url{http://www.jstatsoft.org/v48/i10/}
}
\examples{
# mini example with three Likert items and two subject covariates

# using example data "xmpl" in the package
dsgnmat <- patt.design(xmpl, nitems = 3, resptype="rating",
    ia = TRUE, cov.sel="ALL")
head(dsgnmat)

# fit of Critchlov & Fligner (1991) Salad Dressings Data
pattR.fit(salad, nitems = 4)

# alternatively use glm() with patt.design()
sal <- patt.design(salad, nitems = 4, resptype = "ranking")
glm(y ~ A+B+C+D, data = sal, family = poisson)
}
\keyword{package}
