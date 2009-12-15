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
}
\details{
\tabular{ll}{
Package: \tab prefmod\cr
Type: \tab Package\cr
Version: \tab 0.8-17.2\cr
Date: \tab 2009-12-07\cr
Depends: \tab stats\cr
License: \tab GPL\cr
}

%Index:
%\preformatted{
%cemspc                  Data: CEMS
%checkMIS                check and display NAs in paired comparisons data
%dat4                    Data: dat4
%issp2000                Data: ISSP 2000 Survey on Environmental Issues
%llbt.design             Loglinear Bradley-Terry Model (LLBT) - Design
%                        Matrix Generation
%llbt.fit                Function to fit an LLBT
%baseball                Data: Baseball Games
%music                   Data: Music (Music (General social survey)
%patt.design             Paired Comparison Patterns - Design Matrix
%                        Generation
%pattL.fit               Function to fit a pattern model for ratings
%                        (Likert items)
%pattPC.fit              Function to fit a pattern model for paired
%                        comparisons
%pattR.fit               Function to fit a pattern model for rankings
%salad                   Data: Salad Dressings (Critchlow and Fligner)
%xmpl                    Data: Example Data Set
%}
}
\author{
Reinhold Hatzinger

Maintainer: Reinhold Hatzinger <reinhold.hatzinger@wu-wien.ac.at>
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
