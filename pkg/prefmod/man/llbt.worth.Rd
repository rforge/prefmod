\name{llbt.worth}
\alias{llbt.worth}
\title{Function to calculate and print worth parameters from LLBT model results}
\description{
   Worth parameter are calculated from the results of an LLBT model fit, i.e., from
   \code{\link{llbtPC.fit}}.
   The resulting estimates for all groups are based on the subject covariates as
   specified in the model formula (\code{formel}) of the  LLBT model.
}
\usage{
llbt.worth(obj, obj.names = NULL, outmat = "worth")
}
\arguments{
  \item{obj}{Object of class \code{llbtMod} obtained from an LLBT model fit, i.e.,
        \code{\link{llbtPC.fit}}.}
  \item{obj.names}{names for the objects. If \code{NULL}, the names as specified
        in \code{obj} are used.}
  \item{outmat}{a matrix of estimated worth parameters (\code{outmat = "worth"}, the default) or
        LLBT model parameters (\code{outmat = "lambda"}).}
}
\value{
  \code{llbt.worth} returns a matrix
  of worth or model parameters. If subject covariates have been specified, each column
  represents a groups defined by the crossclassification of the subject covariates.

  The function \code{\link{plotworth}} gives a plot of the estimates.
}
\author{Reinhold Hatzinger}
\seealso{\code{\link{llbtPC.fit}}, \code{\link{plotworth}}}
\examples{
    ## fit only first three objects with SEX effect
    data(cemspc)
    m2<-llbtPC.fit(cemspc, nitems=3, formel=~SEX, elim=~SEX, undec=TRUE)

    ## calculate and print worth parameters
    m2worth<-llbt.worth(m2)
    m2worth
}
\keyword{models}
