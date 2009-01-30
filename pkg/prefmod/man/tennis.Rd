\name{tennis}
\alias{tennis}
\docType{data}
\title{Data (paired comparisons): Preferred Interview Partner}
\description{
  The data describes results from a paired comparison study where 68 male and 96
  female students were asked whom they would prefer to interview. The potential
  interview partners were Bonnie Blair, Jackie Joyner, and Jennifer Capriati.
}
\usage{data(tennis)}
\format{
  A data frame with 16 observations on the following 5 variables.
  \describe{
    \item{\code{n}}{counts of response pattern (\code{C1,C2,C3})}
    \item{\code{C1}}{Blair vs. Joyner (1 \dots Blair preferred, -1 \dots Joyner preferred)}
    \item{\code{C2}}{Blair vs. Capriati (1 \dots Blair preferred, -1 \dots Capriati preferred)}
    \item{\code{C3}}{Joyner vs. Capriati (1 \dots Joyner preferred, -1 \dots Capriati preferred)}
    \item{\code{SEX}}{a numeric vector: 1 \dots male, 2 \dots female}
  }
}
\references{B\"ockenholt , U. and Dillon , W.R., Modelling within-subject
dependencies in ordinal paired comparison data, Psychometrika, 1997, 62, p.411-434
}
\examples{
data(tennis)
tdat<-expand.mat(tennis[,-1],tennis[,1])
head(tdat)
}
\keyword{datasets}
