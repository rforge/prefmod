\name{patt.worth}
\alias{patt.worth}
\title{Function to calculate and print worth parameters from pattern model results}
\description{
   Worth parameter are calculated from the results of a pattern model fit.
   The resulting estimates for all groups are based on the subject covariates as
   specified in the model formula (\code{formel}) of a pattern model.
}
\usage{
patt.worth(obj, obj.names = NULL)
}
\arguments{
  \item{obj}{Object of class \code{pattMod} obtained from pattern model fit.}
  \item{obj.names}{names for the objects, for repeated measurement models just the names of
        objects for the first time point}
}
\value{
  \code{patt.worth} returns an object of class \code{pattW}, which is in fact the matrix .
  of worth parameters for all groups appropriately labelled.
  The function print (i.e., print.pattW) can be used to print the matrix, the function
  \code{\link{plotworth}} gives a plot of the estimates. The matrix
  itself can be obtained by using \code{\link{unclass}}.
}
\author{Reinhold Hatzinger}
\examples{
    ## fit only first three objects with SEX effect
    data(cemspc)
    m2<-pattPC.fit(cemspc, nitems=3, formel=~SEX, elim=~SEX, undec=TRUE)

    ## calculate and print worth parameters
    m2worth<-patt.worth(m2)
    m2worth

    ## display the transpose of the worth parameter matrix
    t(unclass(m2worth))
}
\keyword{models}
