\name{patt.worth}
\alias{patt.worth}
\title{Function to calculate and print worth parameters from pattern model results}
\description{
   Worth parameter are calculated from the results of a pattern model fit, i.e., from
   \code{\link{pattPC.fit}}, \code{\link{pattR.fit}}, \code{\link{pattL.fit}}, and
   \code{\link{pattLrep.fit}}.
   The resulting estimates for all groups are based on the subject covariates as
   specified in the model formula (\code{formel}) of a pattern model.
}
\usage{
patt.worth(obj, obj.names = NULL, outmat = "worth")
}
\arguments{
  \item{obj}{Object of class \code{pattMod} obtained from pattern model fit.}

  \item{obj.names}{names for the objects, for repeated measurement models just the names of
        objects for the first time point}
  \item{outmat}{a matrix of estimated worth parameters (\code{outmat = "worth"}, the default) or
        pattern model parameters (\code{outmat = "lambda"}).}
}
\value{
  \code{patt.worth} returns a matrix
  of worth or model parameters. If subject covariates have been specified, each column
  represents a groups defined by the crossclassification of the subject covariates.

  The function \code{\link{plotworth}} gives a plot of the estimates.
}
\author{Reinhold Hatzinger}
\seealso{\code{\link{pattPC.fit}}, \code{\link{pattR.fit}}, \code{\link{pattL.fit}}, and
   \code{\link{pattLrep.fit}}, \code{\link{plotworth}}}
\examples{
    ## fit only first three objects with SEX effect
    data(cemspc)
    m2<-pattPC.fit(cemspc, nitems=3, formel=~SEX, elim=~SEX, undec=TRUE)

    ## calculate and print worth parameters
    m2worth<-patt.worth(m2)
    m2worth
}
\keyword{models}
