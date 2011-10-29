\name{carconf}
\alias{carconf}
\docType{data}
\title{Data (partial rankings): Car Configurator}
\description{
Online configuration systems allow customers to actively participate in
the creation of products and become increasingly important in various
market sectors. Dabic and Hatzinger report a study on car configurators that aimed
at investigating the effects of certain person characteristics (such as gender)
on the configuration process. Subjects were asked to configure a car according
to their preferences. They could choose freely from several modules: such as
exterior and interior design, technical equipment, brand, price, and producing country.
The order of module choice was recorded as ranks. Since not all modules had to
be chosen the response format was partial rankings.}
\usage{data(carconf)}
\format{
  A data frame with 435 observations on the following 9 variables.
  \describe{
    \item{\code{price}}{rank (1 highest preference)}
    \item{\code{exterior}}{rank}
    \item{\code{brand}}{rank}
    \item{\code{tech.equip}}{rank}
    \item{\code{country}}{rank}
    \item{\code{interior}}{rank}
    \item{\code{sex}}{1 female 2 male}
    \item{\code{age}}{1 17--29 years, 2 30--49 years, 3 50+ years}
    \item{\code{segment}}{preferred car type: 1 premium-class, 2 medium-class, 3 low-budget}
  }
}
\source{
Dabic, M., Hatzinger, R. (2009). Zielgruppenadaequate Ablaeufe in
Konfigurationssystemen ?? eine empirische Studie im Automobilmarkt --
Partial Rankings. In R Hatzinger, R Dittrich, T Salzberger (eds.),
Praeferenzanalyse mit R: Anwendungen aus Marketing, Behavioural
Finance und Human Resource Management. Facultas. }
\examples{
data(carconf)
head(carconf)
}
\keyword{datasets}
