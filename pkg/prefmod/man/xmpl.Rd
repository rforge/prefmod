\name{xmpl}
\alias{xmpl}
\docType{data}
\title{Data (Likert items): Example Data Set}
\description{
  Data to illustrate the usage of \code{\link{patt.design} for Likert items}.
}
\usage{data(xmpl)}
\format{
  A data frame with 100 observations on the following 5 variables.
  \describe{
    \item{\code{I1}}{responses to first Likert item, (1) \emph{strong agreement}, to
                     (5) \emph{strong disagreement}}
    \item{\code{I2}}{responses to second Likert item, like \code{I1}}
    \item{\code{I3}}{responses to third Likert item, like \code{I1}}
    \item{\code{SEX}}{(1) \emph{male}, (2) \emph{female}}
    \item{\code{EDU}}{(1) \emph{low education}, (2) \emph{high education}}
  }
  All values are numeric.
}
\details{
  Datasets and/or dataframes used in \code{\link{patt.design}} are required
  to have the follwing structure:
  \itemize{
     \item All values must be numeric.
     \item The item responses must be in the leftmost columns (such as \code{I1}
           to \code{I3} above).
     \item Categorical subject covariates follow the item responses (rightmost columns)
           and their levels must be specified as consecutive integers. If in a used
           datafile or dataframe these are defined as \pkg{R} factors they will be
           converted to integers. This is not possible if characters are used
           as factor levels and, consequently, \code{patt.design} will produce an error.
  }
}
\examples{
data(xmpl)
des<-patt.design(xmpl, nitems = 3, resptype = "rating",
        cov.sel = "SEX")
head(des)
}
\keyword{datasets}
