\name{patt.worth}
\alias{patt.worth}
\title{Function to calculate and print worth parameters from pattern model results}
\description{
   Worth parameter are calculated from the results of a pattern model fit, i.e., from
   \code{\link{pattPC.fit}}, \code{\link{pattR.fit}}, \code{\link{pattL.fit}}, and
   \code{\link{pattLrep.fit}} or from a \code{gnm}-fit, respectively. For the latter, the
   function only works if the design matrix had been generated using \code{patt.design}.
}
\usage{
patt.worth(fitobj, obj.names = NULL, outmat = "worth")
}
\arguments{
  \item{fitobj}{Object of class \code{pattMod} obtained from pattern model fit.}

  \item{obj.names}{names for the objects, for repeated measurement models just the names of
        objects for the first time point}
  \item{outmat}{a matrix of estimated worth parameters (\code{outmat = "worth"}, the default) or
        pattern model parameters (\code{outmat = "lambda"}).}
}
\details{If the pattern model includes categorical subject covariates, the
  function provides estimates for all groups formed by the full
  crossclassification. Numerical subject covariates are not implemented (yet)(see Warning below).
}
\section{Warning}{
  If the pattern model has been fitted including numerical subject covariates, they are
  ignored. However, estimates for the remaining predictors are calculated for convenience.
  Please note, that these cannot be interpreted as standard estimates but are intercepts
  of the regression model where the objects (or reparameterized objects) are explained
  by one or more numerical subject covariates.
}
\value{
  \code{patt.worth} returns a matrix
  of worth or model parameters. If subject covariates have been specified, each column
  represents a groups defined by the crossclassification of the subject covariates.

  The function \code{\link{plot}} gives a plot of the estimates.
}
\author{Reinhold Hatzinger}
\seealso{\code{\link{pattPC.fit}}, \code{\link{pattR.fit}}, \code{\link{pattL.fit}}, and
   \code{\link{pattLrep.fit}}, \code{\link{plot}}}
\examples{
    ## fit only first three objects with SEX effect
    data(cemspc)
    m2<-pattPC.fit(cemspc, nitems=3, formel=~SEX, elim=~SEX, undec=TRUE)

    ## calculate and print worth parameters
    m2worth<-patt.worth(m2)
    m2worth
}
\keyword{models}
