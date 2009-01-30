\name{cemspc}
\alias{cemspc}
\docType{data}
\title{Data (paired comparisons with undecided): CEMS (Community of European management schools)}
\description{
A survey of 303 students was carried out to examine the student's
preferences of 6 universities (London, Paris, Milano, St.Gallen,
Barcelona and Stockholm). The first 15 variables of the 17 indicate
the preferences of the subjects. For a given comparison the reponses
are coded by 1 if the first university was preferred, 2 if the second
university was preferred and 0 if no decision was made. The variable
ENG characterises the knowledge of English and the variable
SEX characterises the gender.
}
\usage{data(cemspc)}
\format{
  A data frame with 303 observations on the following 17 variables.
  \describe{
    \item{\code{V1}}{comparison of London to Paris}
    \item{\code{V2}}{comparison of London to Milano}
    \item{\code{V3}}{comparison of Paris to Milano}
    \item{\code{V4}}{comparison of London to St.Gallen}
    \item{\code{V5}}{comparison of Paris to St.Gallen}
    \item{\code{V6}}{comparison of Milano to St.Gallen}
    \item{\code{V7}}{comparison of London to Barcelona}
    \item{\code{V8}}{comparison of Paris to Barcelona}
    \item{\code{V9}}{comparison of Milano to Barcelona}
    \item{\code{V10}}{comparison of St.Gallen to Barcelona}
    \item{\code{V11}}{comparison of London to Stockholm}
    \item{\code{V12}}{comparison of Paris to Stockholm}
    \item{\code{V13}}{comparison of Milano to Stockholm}
    \item{\code{V14}}{comparison of St.Gallen to Stockholm}
    \item{\code{V15}}{comparison of Barcelona to Stockholm}
    \item{\code{ENG}}{Knowledge of English : 1= good, 2= poor}
    \item{\code{SEX}}{Gender : 1= female, 2= male}
  }
}
\references{
  R. Dittrich, R.Hatzinger and W.Katzenbeisser: Modelling the effect of subject-specific covariates in paired comparison studies with an application to university rankings. Applied Statistics, Vol. 47 (1998), Part 4
}
\examples{
data(cemspc)
par(mfrow=c(4,4))
for (i in 1:15) barplot(table(cemspc[,i]))
par(mfrow=c(1,1))
}
\keyword{datasets}
