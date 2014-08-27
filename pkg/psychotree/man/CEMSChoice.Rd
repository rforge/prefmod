\name{CEMSChoice}
\alias{CEMSChoice}
\encoding{latin1}

\title{CEMS University Choice Data}

\description{
  Preferences of 303 students from WU Wien for different CEMS universities.
}

\usage{data("CEMSChoice")}

\format{
A data frame containing 303 observations on 10 variables.
\describe{
  \item{preference}{Paired comparison of class \code{\link{paircomp}}.
    Preferences for all 15 paired comparisons from 6 objects:
    London, Paris, Milano, St. Gallen, Barcelona, Stockholm.}
  \item{study}{Factor coding main discipline of study: commerce,
    or other (economics, business administration, business education).}
  \item{english}{Factor coding knowledge of English (good, poor).}
  \item{french}{Factor coding knowledge of French (good, poor).}
  \item{spanish}{Factor coding knowledge of Spanish (good, poor).}
  \item{italian}{Factor coding knowledge of Italian (good, poor).}
  \item{work}{Factor. Was the student working full-time while studying?}
  \item{gender}{Factor coding gender.}
  \item{intdegree}{Factor. Does the student intend to take an international degree?}
  \item{preference1998}{Paired comparison of class \code{\link{paircomp}}.
    This is like \code{preference} but the comparisons between Barcelona
    an Stockholm are (erroneously) reversed, see below.}
}
}

\details{
  Students at Wirtschaftsuniversität Wien (\url{http://www.wu.ac.at/}) can study abroad
  visiting one of currently 17 CEMS universities (Community of European Management
  Schools and International Companies). Dittrich et al. (1998) conduct and analyze
  a survey of 303 students to examine the student's preferences for 6 universities:
  London School of Economics, HEC Paris, Università Commerciale Luigi Bocconi (Milano),
  Universität St. Gallen, ESADE (Barcelona), Handelshögskolan i Stockholm.
  To identify reasons for the preferences, several subject covariates (including
  foreign language competence, gender, etc.) have been assessed. Furthermore,
  several object covariates are attached to \code{preference} (and \code{preference1998}):
  the universities' field of \code{specialization} (economics, management science, finance)
  and location (Latin country, or other).

  The correct data are available in the online complements to Dittrich et al. (1998).
  However, the accompanying analysis was based on an erroneous version of the data in which the
  choices for the last comparison pair (Barcelona : Stockholm) were accidentally reversed. See the
  corrigendum in Dittrich et al. (2001) for further details. The variable \code{preference}
  provides the correct data and can thus be used to replicate the analysis from the
  corrigendum (Dittrich et al. 2001). For convenience, the erroneous version is
  provided in \code{preference1998} which can therefore be used to replicate the
  (incorrect) original analysis (Dittrich et al. 1998).
}

\source{
The Royal Statistical Society Datasets Website.

\url{http://www.blackwellpublishing.com/rss/Readmefiles/dittrich.htm}
}

\references{
Dittrich R, Hatzinger R, Katzenbeisser W (1998).
  Modelling the Effect of Subject-Specific Covariates in Paired Comparison Studies
  with an Application to University Rankings,
  \emph{Journal of the Royal Statistical Society C}, \bold{47}, 511--525.

Dittrich R, Hatzinger R, Katzenbeisser W (2001).
  Corrigendum: Modelling the Effect of Subject-Specific Covariates in Paired Comparison
  Studies with an Application to University Rankings,
  \emph{Journal of the Royal Statistical Society C}, \bold{50}, 247--249.
}

\seealso{\code{\link{paircomp}}}

\examples{
data("CEMSChoice", package = "psychotree")
summary(CEMSChoice$preference)
covariates(CEMSChoice$preference)
}

\keyword{datasets}
