\name{music}
\alias{music}
\docType{data}
\title{Data (ratings): Music (General social survey)}
\description{
  The General Social Survey has been conducted by the National Opinion Research Center in 1993. The dataset focusses on 11 items with a 5-point Likert type response scale. The respondents were asked about their music preferences.
}
\usage{data(music)}
\format{
  A data frame with 1500 observations on the following 11 variables. The variables are items with the question for liking bigband/ blugrass,.. very much, like it, have mixed feelings, dislike it and dislike it very much. Missing values are categorized as NA, so it is not a numeric data frame.
  \describe{
    \item{\code{bigband}}{responses to bigband item}
    \item{\code{blugrass}}{responses to blugrass item}
    \item{\code{country}}{responses to country item}
    \item{\code{blues}}{responses to blues item}
    \item{\code{musicals}}{responses to musicals item}
    \item{\code{klassik}}{responses to klassik item}
    \item{\code{folk}}{responses to folk item}
    \item{\code{jazz}}{responses to jazz item}
    \item{\code{opern}}{responses to opern item}
    \item{\code{rap}}{responses to rap item}
    \item{\code{hvymetal}}{responses to hvymetal item}
  }
}
\details{
  The 11 variables are items to be anserwed on a 5-point Likert type scale with response categories: (1) like it very much, (2) like it, (3) have mixed feelings, (4) dislike it and (5) dislike it very much.
}
\source{
    Dataset: General Social Survey 1993 [United States]\cr
    CPANDA Identification Number: a00006\cr
    Cultural Policy and the Arts National Data Archive\cr
   \url{http://www.cpanda.org/cpanda/getDDI.xq?studyID=a00006#dataDscr}\cr
}
\references{
This is just an abstract for further information concerning the General Social Survey dataset have a look at the Cultural Policy and the Arts National Data Archive.
}
\examples{
data(music)
summary(music)
}
\keyword{datasets}
