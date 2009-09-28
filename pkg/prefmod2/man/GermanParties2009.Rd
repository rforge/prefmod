\name{GermanParties2009}
\alias{GermanParties2009}
\encoding{latin1}

\title{Choice among German Political Parties}

\description{
  Preferences of 192 respondents choosing among five German political parties
  and abstention from voting.
}

\usage{data("GermanParties2009")}

\format{
A data frame containing 192 observations on 6 variables.
\describe{
  \item{preference}{Paired comparison of class \code{\link{paircomp}}.
    All 15 pairwise choices among five German parties and abstention from
    voting.}
  \item{ordered.pref}{Ordered paired comparison of class
    \code{\link{paircomp}}. Same as \code{preference}, but within-pair order
    is recognized.}
  \item{gender}{Factor coding gender.}
  \item{age}{Integer. Age of the respondents in years.}
  \item{eduction}{Ordered factor. Level of education: 1 no degree, 2
    Hauptschule (Secondary General School), 3 Realschule (Intermediate
    Secondary School), 4 Gymnasium (High School), 5 Studium (University)}
  \item{crisis}{Factor. Do you feel affected by the economic crisis?}
}
}

\details{
  A survey was conducted at the Department of Psychology, Universität
  Tübingen, in June 2009, three months before the German election. The sample
  was stratified by gender and age (younger versus older than 30 years) with
  48 participants in each group.

  The parties to be compared were Die Linke (socialists), Die Grünen
  (ecologists), SPD (social democrats), CDU/CSU (conservatives), and FDP
  (liberals). In addition, there was the option of abstaining from voting
  (coded as \code{none}).

  Participants were presented with 15 pairs of options in random order. On
  each trial, their task was to choose the party they would rather vote for at
  an election for the German parliament. A given participant compared each
  pair in only one order, hence the NA's in \code{ordered.pref}.

  In order to minimize response biases, the pairs of options were read to the
  participants one at a time. Participants made their choices by crossing
  either \dQuote{First Option} or \dQuote{Second Option} on an anonymous
  response sheet.
}

% \references{
%   FIX ME
% }

\seealso{\code{\link{paircomp}}}

\examples{
data("GermanParties2009", package = "prefmod2")
summary(GermanParties2009$preference)
}

\keyword{datasets}
