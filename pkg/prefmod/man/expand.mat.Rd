\name{expand.mat}
\alias{expand.mat}
\title{Utility function to expand aggregated data}
\description{
  The function expands aggregated data into casewise data. For instance, for a
  contingency table given in the form of a design matrix and corresponding counts
  the function sets up a matrix where each design row is repeated according to
  the frequencies for that row.
}
\usage{
expand.mat(mat, freq)
}
\arguments{
  \item{mat}{a matrix (or column vector) or data frame to be expanded}
  \item{freq}{a vector of counts}
}
%\details{}
\value{
  the expanded matrix
}
\note{This utility allows to generate input data for the for the design generating
      and model fitting functions of the \code{prefmode} package
      from aggregated data.
}
\examples{
data(tennis)
tdat<-expand.mat(tennis[,-1],tennis[,1])
head(tdat)
}
\keyword{utilities}
