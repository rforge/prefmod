\name{SPISA}
\alias{SPISA}

\title{SPIEGEL Studentenpisa Data (Subsample)}

\description{
  A subsample from the general knowledge quiz \dQuote{Studentenpisa}
  conducted online by the German weekly news magazine SPIEGEL. The
  data contain the quiz results from 45 questions as well as sociodemographic
  data for 1075 university students from Bavaria.
}

\usage{data("SPISA")}

\format{
A data frame containing 1075 observations on 6 variables.
\describe{
  \item{spisa}{matrix with \code{0}/\code{1} results from
    45 questions in the quiz (indicating wrong/correct answers).}
  \item{gender}{factor indicating gender.}
  \item{age}{age in years.}
  \item{semester}{numeric indicating semester of university enrollment.}
  \item{elite}{factor indicating whether the university the student is
    enrolled in has been granted \dQuote{elite} status by the German
    \dQuote{excellence initiative}.}
  \item{spon}{ordered factor indicating frequency of accessing the
    SPIEGEL online (SPON) magazine.}
}
}

\details{
  An online quiz for testing one's general knowledge was conducted by
  the German weekly news magazine SPIEGEL in 2009. Overall, about 700,000
  participants answered the quiz and a set of sociodemographic questions.
  The general knowledge quiz consisted of a total of 45 items from five
  different topics: politics, history, economy, culture and natural sciences.
  For each topic, four different sets of nine items were available, that
  were randomly assigned to the participants. A thorough analysis and
  discussion of the original data set is provided in Trepte and Verbeet (2010).
  
  Here, we provide the subsample of university students enrolled in the
  federal state of Bavaria, who had been assigned questionnaire number 20
  (so that all subjects have answered the same set of items). Excluding
  all incomplete records, this subsample contains 1075 observations.
  
  The data are analyzed in Strobl et al. (2010) whose analysis is replicated
  in \code{vignette("raschtree", package = "psychotree")}.
}

\references{
  Strobl, C., Kopf, J., and Zeileis, A. (2010). A New Method for Detecting
  Differential Item Functioning in the Rasch Model. Unpublished Manuscript.

  SPIEGEL Online (2009). Studentenpisa -- Alle Fragen, alle Antworten.
  In German. Accessed 2010-10-26.
  \url{http://www.spiegel.de/unispiegel/studium/0,1518,620101,00.html}
  
  Trepte, S., and Verbeet, M. (2010). Allgemeinbildung in
  Deutschland -- Erkenntnisse aus dem SPIEGEL-Studentenpisa-Test.
  ISBN 978-3-531-17218-7. VS Verlag, Wiesbaden.
}

\seealso{\code{\link{raschtree}}}

\examples{
## data
data("SPISA", package = "psychotree")

## summary of covariates
summary(SPISA[,-1])

## histogram of raw scores
hist(rowSums(SPISA$spisa), breaks = 0:45 + 0.5)

\dontrun{
## See the following vignette for a tree-based DIF analysis
vignette("raschtree", package = "psychotree")
}
}

\keyword{datasets}
