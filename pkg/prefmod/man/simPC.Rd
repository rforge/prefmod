\name{simPC}
\alias{simPC}
\alias{simR}
\title{Utility function to simulate paired comparison or rankings data
}
\description{The function generates a random paired comparison
            data matrix (two response categories, no undecided) or a rankings data matrix
            optionally based on user specified worth parameters.
}
\usage{
simPC(nobj, nobs, worth = NULL, seed = NULL, pr = FALSE)
simR(nobj, nobs, worth = NULL, seed = NULL, pr = FALSE)
}
\arguments{
  \item{nobj}{Number of objects.
}
  \item{nobs}{Number of cases.
}
  \item{worth}{
     If \code{NULL}, values are drawn from uniform distribution (interval from 0 to 1), otherwise the user may specify
     arbitrary values. In both cases, the values are transformed into [0,1] and scaled to sum up to 1.
}
  \item{seed}{Starting value for the random number generator.
}
  \item{pr}{If \code{TRUE}, the used worth parameters are printed.
}
}
\value{The random data matrix as a data frame
}
\author{Reinhold Hatzinger
}
\examples{
data<-simPC(4,10,worth=1:4,seed=123456)
data
}
\keyword{datagen}
