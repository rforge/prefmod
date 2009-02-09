\name{TopModel2007}
\alias{TopModel2007}
\encoding{latin1}

\title{Attractiveness of Germany's Next Topmodels}

\description{
  Preferences of 192 respondents judging the attractiveness of the top six
  contestants of the TV show Germany's Next Topmodel 2007.
}

\usage{data("TopModel2007")}

\format{
A data frame containing 192 observations on 6 variables.
\describe{
  \item{preference}{Paired comparison of class \code{\link{paircomp}}.
    Preferences for all 15 paired comparisons from 6 contestants:
    Barbara, Anni, Hana, Fiona, Mandy, and Anja.}
  \item{gender}{Factor coding gender.}
  \item{age}{Integer. Age of the respondents in years.}
  \item{q1}{Factor. Do you recognize the women on the pictures?/Do you know
    the TV show Germany's Next Topmodel?}
  \item{q2}{Factor. Did you watch Germany's Next Topmodel regularly?}
  \item{q3}{Factor. Did you watch the final show of Germany's Next
    Topmodel?/Do you know who won Germany's Next Topmodel?}
}
}

\details{
  The survey was conducted in 2007 shortly after the final show. The sample
  was stratified by gender and age (younger versus older than 30 years) with
  48 participants in each group.

  Digital photographs (resolution 303 times 404 pixels) of the top six
  contestants were available from the ProSieben web page
  (\url{http://www.prosieben.de/}) at the time of the survey. The photos were
  selected to be comparable, showing the contestant's face and the upper part
  of the body, all women being casually dressed.

  Participants were presented with all 15 pairs of photographs. On each trial,
  their task was to judge which of the two women on the photos was the more
  attractive. In order to assess the participants' expertise, additional
  questions regarding their familiarity with the show were asked after the
  pairwise comparisons were completed.

  The actual ranking, as resulting from sequential elimination during the
  course of the show, was (from first to last place): Barbara, Anni, Hana,
  Fiona, Mandy, Anja.
}

\references{
Germany's Next Topmodel. In Wikipedia, The Free Encyclopedia.

\url{http://en.wikipedia.org/w/index.php?title=Germany%27s_Next_Topmodel&oldid=268802871}
}

\seealso{\code{\link{paircomp}}}

\examples{
data("TopModel2007", package = "prefmod2")
summary(TopModel2007$preference)
xtabs(~ gender + I(age < 30), TopModel2007)
}

\keyword{datasets}
