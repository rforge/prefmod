\name{trdel}
\alias{trdel}
\docType{data}
\title{Data (paired comparisons): Training delivery modes}
\description{
The dataset \code{trdel} contains data from a paired comparison
study to investigate which of five training delivery modes trainees prefer
(Schoell and Veith, 2011).
The modes were computer-based (\code{CO}), TV-based (\code{TV}), paper-based (\code{PA}),
audio-based (\code{AU}) and classroom-based (\code{CL}) training. Study
participants were unemployed persons in the labour market training of
the Austrian labour market service (AMS). To account for trainee
characteristics that might affect the preference order the variables
gender, age, and learning personality type were recorded. These variables were
coded as \code{sex} (1 male, 2 female), \code{age} (numeric in years), \code{ltype} (1
accomodator, 2 diverger, 3 converger, 4 assimilator). The learning personality
types were identified from a questionnaire.
}
\usage{data(trdel)}
\format{
  A data frame with 198 observations on the following 14 variables.
  \describe{
    \item{\code{V1,V2,V3,V4,V5,V6,V7,V8,V9,V10}}{paired comparisons in standard order: \code{CO:TV}, \code{CO:PA}, etc.
                           1 first object preferred, 2 second object preferred}
    \item{\code{ltype}}{learning types: 1 accomodator, 2 diverger, 3 converger, 4 assimilator}
    \item{\code{age}}{numeric in years}
    \item{\code{sex}}{1 male, 2 female}
  }
}
\source{
Schoell, B., Veith, S. (2011). Learning style evaluation and
preferred training delivery modes in labour market training (in
German). Master's thesis, Vienna University of Economics and
Business. }
\examples{
data(trdel)
head(trdel)
}
\keyword{datasets}
