\name{baseball}
\alias{baseball}
\docType{data}
\title{Data (paired comparisons): Baseball Games}
\description{
  The result of the 1987 season for seven baseball teams in the Eastern Division of the American League according to the (home team, away team) classification are shown.
}
\usage{data(baseball)}
\format{
  Baseball is a numeric vector with the results for the season according to the (home team, away team) classification.
 }
\details{
The results of the seven teams Milwaukee, Detroit, Toronto, New York, Boston, Cleveland and Baltimore, that play 13 games each. There is no possibility of ending in a draw.
}
\references{
   Alan Agresti, \emph{Categorical Data Analysis (Second Edition)}, 2002 pages 437 and 438

   R. Dittrich, R. Hatzinger, and W. Katzenbeisser, Fitting paired comparison models in GLIM.
   \emph{GLIM newsletter} 1997

   }
\examples{
data(baseball)
str(baseball)
}
\keyword{datasets}
