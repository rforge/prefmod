\encoding{UTF-8}
\name{xmpl}
\alias{xmpl}
\docType{data}
\title{Data (Likert items): Example Data Set}
\description{
  Data to illustrate the usage of \code{\link{patt.design} for rating scale (Likert type) items}.
}
\usage{xmpl}
\format{%
  A data frame with 100 observations on 5 numeric variables.
  The first three variables (\code{I1}, \code{I2}, \code{I3}) are the rating scale (Likert type) items with 5 response categories, ranging from 1 (strong agreement) to 5 (strong disagreement).

  \describe{%
    \item{\code{I1}}{response to item 1}
    \item{\code{I2}}{response to item 2}
    \item{\code{I3}}{response to item 3}
    \item{\code{SEX}}{(1) \emph{male}, (2) \emph{female}}
    \item{\code{EDU}}{(1) \emph{low education}, (2) \emph{high education}}
  }
}
\details{%
  Datasets in data files or Data frames used in \code{\link{patt.design}} require the following structure:
  \itemize{%
    \item All values must be numeric.
    \item The item responses must be in the leftmost columns (such as \code{I1} to \code{I3} above).
    \item Categorical subject covariates follow the item responses (rightmost columns) and their levels must be specified as consecutive integers.
      If in a used datafile or dataframe these are defined as \pkg{R} factors they will be converted to integers.
      This is not possible if characters are used as factor levels and, consequently, \code{patt.design} will produce an error.
  }
}
\examples{
des <- patt.design(xmpl, nitems = 3, resptype = "rating", cov.sel = "SEX")
head(des)
}
\keyword{datasets}
