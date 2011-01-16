\name{cemspc}
\alias{cemspc}
\docType{data}
\title{Data (paired comparisons with undecided): CEMS (Community of European management schools)}
\description{
A survey of 303 students was carried out to examine the student's
preferences of 6 universities (London, Paris, Milano, St.Gallen,
Barcelona and Stockholm) with a 17 items questionnaire. The first 15 variables indicate
the subjects' preferences. For a given comparison the responses
are coded as 0 if the first university was preferred, 2 if the second
university was preferred and 1 if no decision was made. The variable
ENG contains the knowledge of English and the variable
SEX contains the gender.
}
\usage{data(cemspc)}
\format{
  A data frame with 303 rows for the subjects containing the outcome of the 15 comparisons and the two covariates 
  \describe{
    \item{\code{V1}}{London vs. Paris}
    \item{\code{V2}}{London vs. Milano}
    \item{\code{V3}}{Paris vs. Milano}
    \item{\code{V4}}{London vs. St.Gallen}
    \item{\code{V5}}{Paris vs. St.Gallen}
    \item{\code{V6}}{Milano vs. St.Gallen}
    \item{\code{V7}}{London vs. Barcelona}
    \item{\code{V8}}{Paris vs. Barcelona}
    \item{\code{V9}}{Milano vs. Barcelona}
    \item{\code{V10}}{St.Gallen vs. Barcelona}
    \item{\code{V11}}{London vs. Stockholm}
    \item{\code{V12}}{Paris vs. Stockholm}
    \item{\code{V13}}{Milano vs. Stockholm}
    \item{\code{V14}}{St.Gallen vs. Stockholm}
    \item{\code{V15}}{Barcelona vs. Stockholm}
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
