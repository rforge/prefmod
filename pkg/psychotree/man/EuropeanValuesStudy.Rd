\name{EuropeanValuesStudy}
\alias{EuropeanValuesStudy}

\title{European Values Study}

\description{
  A sample of the 1999 European Values Study (EVS) containing an assessment of
  materialism/postmaterialism in 3584 respondents from 32 countries.
}

\usage{data("EuropeanValuesStudy")}

\format{
A data frame containing 3584 observations on 10 variables.
\describe{
  \item{country}{Factor coding the country of a respondent.}
  \item{gender}{Factor coding gender.}
  \item{birthyear}{Numeric. Year of birth.}
  \item{eduage}{Numeric. Age when full time education was or will be
    completed.}
  \item{marital}{Factor. Current legal marital status.}
  \item{employment}{Ordered factor. Employment and number of working hours.}
  \item{occupation}{Factor. What is/was your main job?}
  \item{income}{Ordered factor. Income of household in ten categories from
    10 percent lowest to 10 percent highest income category.}
  \item{paircomp}{Paired comparison of class \code{\link{paircomp}}.
    Five pairwise choices among four important political goals derived from
    a double-choice task (see Details).}
  \item{county2}{Factor. Country group according to postmaterialism
    (see Details).}
}
}

\details{
  The data are part of a larger survey conducted in 1999 in 32 countries in
  Europe (see \url{http://www.europeanvaluesstudy.eu/}). Vermunt (2003)
  obtained a sample from 10 percent of the available cases per country,
  yielding 3584 valid cases.

  The item in the 1999 European Values Study questionnaire aiming at recording
  materialism/postmaterialism reads as follows:

  There is a lot of talk these days about what the aims of this country
  should be for the next ten years. On this card are listed some of the goals
  which different people would give top priority. If you had to choose, which
  of the things on this card would you say is most important? And which would
  be the next most important?
  
  A Maintaining order in the nation\cr
  B Giving people more say in important government decisions\cr
  C Fighting rising prices\cr
  D Protecting freedom of speech

  The double-choice task implies a partial ranking of the alternatives and
  (assuming transitivity) an incomplete set of paired comparisons for each
  respondent.

  The country group according to postmaterialism was derived by Vermunt (2003)
  using a latent class model, and applied by Lee and Lee (2010) in a tree
  model.
}

\source{
  Latent GOLD Sample Data Sets Website.

  \url{http://www.statisticalinnovations.com/products/choice_datasets.html}
}

\references{
  Lee, P.H., & Yu, P.L.H. (2010). Distance-based tree models for ranking
  data. \emph{Computational Statistics and Data Analysis}, \bold{54},
  1672--1682.

  Vermunt, J.K. (2003). Multilevel latent class models.
  \emph{Sociological Methodology}, \bold{33}, 213--239.
}

\seealso{\code{\link{paircomp}}}

\examples{
##
data("EuropeanValuesStudy", package = "psychotree")
summary(EuropeanValuesStudy$paircomp)

## Bradley-Terry tree resulting in similar results compared to
## the (different) tree approach of Lee and Lee (2010)
evs <- na.omit(EuropeanValueStudy)
bt <- bttree(paircomp ~ gender + age + birthyear + marital + employment + income + country2,
  data = evs, alpha = 0.01)
plot(bt, abbreviate = 2)
}

\keyword{datasets}
