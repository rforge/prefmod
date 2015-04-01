\name{SPISA}
\alias{SPISA}
\encoding{latin1}

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
  
  The data are analyzed in Strobl et al. (2010), whose analysis is replicated
  in \code{vignette("raschtree", package = "psychotree")}.

  The full list of items in questionnaire 20 is given below.
  
  Politics:\cr
  Who determines the rules of action in German politics according to the constitution?
  -- The Bundeskanzler (federal chancellor).\cr
  What is the function of the second vote in the elections to the German Bundestag (federal parliament)?
  -- It determines the allocation of seats in the Bundestag.\cr
  How many people were killed by the RAF (Red Army Faction)?
  -- 33.\cr
  Where is Hessen (i.e., the German federal country Hesse) located?
  -- (Indicate location on a map.)\cr
  What is the capital of Rheinland-Pfalz (i.e., the German federal country Rhineland-Palatinate)?
  -- Mainz.\cr
  Who is this? 
  -- (Picture of Horst Seehofer.)\cr
  Which EU institution is elected in 2009 by the citizens of EU member countries?
  -- European Parliament.\cr
  How many votes does China have in the UNO general assembly?
  -- 1.\cr
  Where is Somalia located?
  -- (Indicate location on a map.)

  History:\cr
  The Roman naval supremacy was established through...
  -- ... the abolition of Carthage.\cr
  In which century did the Thirty Years' War take place?
  -- The 17th century.\cr
  Which form of government is associated with the French King Louis XIV?
  -- Absolutism.\cr
  What island did Napoleon die on in exile?
  -- St. Helena.\cr
  How many percent of the votes did the NSDAP receive in the 1928 elections of the German Reichstag?
  -- About 3 percent.\cr
  How many Jews were killed by the Nazis during the Holocaust?
  -- About 6 Million.\cr
  Who is this?
  -- (Picture of Johannes Rau, former German federal president.)\cr
  Which of the following countries is not a member of the EU?
  -- Croatia.\cr
  How did Mao Zedong expand his power in China?
  -- The Long March.

  Economy:\cr
  Who is this?
  -- (Picture of Dieter Zetsche, CEO of Mercedes-Benz.)\cr
  What is the current full Hartz IV standard rate (part of the social welfare) for adults?
  -- 351 Euro.\cr
  What was the average per capita gross national product in Germany in 2007?
  -- About 29,400 Euro.\\
  What is a CEO?
  -- A Chief Executive Officer.\cr
  What is the meaning of the hexagonal ``organic'' logo?
  -- Synthetic pesticides are prohibited.\cr
  Which company does this logo represent?
  -- Deutsche Bank.\cr
  Which German company took over the British automobile manufacturers Rolls-Royce?
  -- BMW.\cr
  Which internet company took over the media group Time Warner?
  -- AOL.\cr
  What is the historic meaning of manufacturies?
  -- Manufacturies were the precursors of industrial mass production.\cr

  Culture:\cr
  Which painter created this painting?
  -- Andy Warhol.\cr
  What do these four buildings have in common?
  -- All four were designed by the same architects.\cr
  Roman numbers: What is the meaning of CLVI?
  -- 156.\cr
  What was the German movie with the most viewers since 1990?
  -- Der Schuh des Manitu.\cr
  In which TV series was the US president portrayed by an African American actor for a long time?
  -- 24.\cr
  What is the name of the bestselling novel by Daniel Kehlmann?
  -- Die Vermessung der Welt (Measuring The World).\cr
  Which city is the setting for the novel \sQuote{Buddenbrooks}?
  -- Lübeck.\cr
  In which city is this building located?
  -- Paris.\cr
  Which one of the following operas is not by Mozart?
  -- Aida.

  Natural sciences:\cr
  Why does an ice floe not sink in the water?
  -- Due to the lower density of ice.\cr
  What is ultrasound not used for?
  -- Radio.\cr
  Which sensory cells in the human eye make color vision possible?
  -- Cones.\cr
  What is also termed Trisomy 21?
  -- Down syndrome.\cr
  Which element is the most common in the Earth's atmosphere?
  -- Nitrogen.\cr
  Which kind of tree does this leaf belong to?
  -- Maple.\cr
  Which kind of bird is this?
  -- Blackbird.\cr
  Where is the stomach located?
  -- (Indicate location on a map of the body.)\cr
  What is the sum of interior angles in a triangle?
  -- 180 degrees.
}

\references{
  Strobl C, Kopf J, Zeileis A (2015). Rasch Trees:
  A New Method for Detecting Differential Item Functioning in the Rasch
  Model. \emph{Psychometrika}, forthcoming. doi:10.1007/s11336-013-9388-3

  SPIEGEL Online (2009). Studentenpisa -- Alle Fragen, alle Antworten.
  In German. Accessed 2010-10-26.
  \url{http://www.spiegel.de/unispiegel/studium/0,1518,620101,00.html}
  
  Trepte S, Verbeet M (2010). Allgemeinbildung in
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
