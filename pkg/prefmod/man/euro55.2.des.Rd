\encoding{UTF-8}
\name{euro55.2.des}
\Rdversion{1.1}
\alias{euro55.2.des}
\docType{data}
\title{
Design data frame for a paired comparison pattern model for rankings
(Eurobarometer 55.2)
}
\description{
Eurobarometer public opinion surveys have been carried out in all
member states of the European Union since 1973. Eurobarometer 55.2
was a special survey collected in 2001 and designed to elicit information
on European experience and perception of science and technology. One
question, of identical form in all states, asked respondents about
their sources of information about science, and requested them to
rank them in order of importance. The items were:
Television, Radio, Press, Scientific Magazines, Internet, School and
University.
}
\usage{euro55.2.des}
\format{
The data is a design data frame generated from the original data using
\code{\link{patt.design}}. Each row represents a (derived) paired comparison
response pattern crossclassified by the categorical subject covariates
\code{SEX} and \code{AGE4}. The columns are:
 \describe{
  \item{\code{y}     }{counts}
  \item{\code{TV}    }{design vector for Television}
  \item{\code{RAD}   }{design vector for Radio}
  \item{\code{NEWSP} }{design vector for Newspapers and magazine}
  \item{\code{SCIMAG}}{design vector for Scientific magazines}
  \item{\code{WWW}   }{design vector for The internet}
  \item{\code{EDINST}}{design vector for School/University}
  \item{\code{SEX}   }{Gender - Factor with 2 levels: (1) male, (2) female}
  \item{\code{AGE4}  }{Age - Factor with 4 levels: (1) 15-24, (2) 25-39, (3) 40-54, (4) 55+}
 }
}
\details{The original data contained  16 029 cases, after removal of cases with missing values
   and improper rankings this design data frame is based on 12216 observations.}
\source{\url{http://ec.europa.eu/research/press/2001/pr0612en.html}
}
\references{Christensen, T. (2001). Eurobarometer 55.2: Europeans, science and technology. Technical
report, European Opinion Research Group, Commission of the European Communities,
Brussels.
%Christensen, Thomas. EUROBAROMETER 55.2: SCIENCE AND TECHNOLOGY, AGRICULTURE, THE EURO, AND INTERNET ACCESS, MAY-JUNE 2001 [Computer file].
%ICPSR03341-v2. Brussels: European Opinion Research Group EEIG [producer], 2001.
%Cologne, Germany: Zentralarchiv fur Empirische Sozialforschung/Ann Arbor,
%MI: Inter-university Consortium for Political and Social Research [distributors], 2006-08-03.
%doi:10.3886/ICPSR03341
}
\examples{str(euro55.2.des)}
\keyword{datasets}
