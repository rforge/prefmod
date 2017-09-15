\encoding{UTF-8}
\name{pattLrep.fit}
\alias{pattLrep.fit}
\title{Function to fit a pattern model for repeated ratings (Likert items)}
\description{%
  Function to fit a pattern model for repeated ratings/Likert items (transformed to paired comparisons) allowing for missing values using a CL approach.
}
\usage{%
pattLrep.fit(obj, nitems, tpoints = 1, formel = ~1, elim = ~1,
        resptype = "ratingT",  obj.names = NULL, undec = TRUE, ia = FALSE,
        iaT = FALSE, NItest = FALSE, pr.it = FALSE)
}
\arguments{
  \item{obj}{either a dataframe or the path/name of the datafile to be read.}
  \item{nitems}{the number of items at one time point.}
  \item{tpoints}{the number of time points.}
  \item{formel}{the formula for subject covariates to fit different preference scales for the objects (see below).}
  \item{elim}{%
    the formula for the subject covariates that specify the table to be analysed.
    If omitted and \code{formel} is not \code{~1} then \code{elim} will be set to the highest interaction between all terms contained in \code{formel}.
    If \code{elim} is specified, the terms must be separated by the \code{*} operator.
  }
  \item{resptype}{is \code{"ratingT"} by default and is reserved for future usage. Any other specification will not change the behaviour of \code{pattLrep.fit}}
  \item{obj.names}{character vector with names for objects.}
  \item{undec}{for paired comparisons with a undecided/neutral category, a common parameter will be estimated if \code{undec = TRUE}.}
  \item{ia}{for each time point interaction parameters between comparisons that have one object in common if \code{ia = TRUE}.}
  \item{iaT}{if \code{iaT = TRUE}, dependence parameters for each item between two successive time points.}
  \item{NItest}{separate estimation of object parameters for complete and incomplete patterns if \code{NItest = TRUE}. Currently, \code{NItest} is set to \code{FALSE} if subject covariates are specified.}
  \item{pr.it}{a dot is printed at each iteration cycle if set to \code{TRUE}}
}
\details{%
  Models including categorical subject covariates can be fitted using the \code{formel} and \code{elim} arguments.
  \code{formel} specifies the actual model to be fitted.
  For instance, if specified as \code{formel = ~SEX} different preference scale for the objects will be estimated for males and females.
  For two or more covariates, the operators \code{+} or \code{*} can be used to model main or interaction effects, respectively.
  The operator \code{:} is not allowed (redundant terms are removed automatically).
  See also \code{\link{formula}}.

  The specification for \code{elim} follows the same rules as for \code{formel}.
  However, \code{elim} specifies the basic contingency table to be set up but does not specify any covariates to be fitted.
  This is done using \code{formel}.
  If, e.g., \code{elim = ~SEX} but \code{formel = ~1}, then the table is set up as if \code{SEX} would be fitted but only one global preference scale is computed.
  This feature allows for the successive fitting of nested models to enable the use of deviance differences for model selection (see example below).
}
\section{Input Data}{
  The input data must have the following order (from left to right): all items at first time point, all items at second time point (with the same order as before), etc. for the other time points, optional subject covariates.
  The responses have to be coded as consecutive integers starting with 1 (or 0).
  The value of 1 (0) means highest \sQuote{endorsement} (agreement) according to the underlying scale.
  Missing values are coded as \code{NA}, rows with less than 2 valid responses are removed from the fit and a message is printed.

  Optional subject covariates have to be specified such that the categories are represented by consecutive integers starting with 1.
  Rows with missing values for subject covariates are removed from the data and a message is printed.
  Again, the leftmost columns in the data must be the ratings, optionally followed by columns for categorical subject covariates.

  The data specified via \code{obj} are supplied using either a data frame or a datafile in which case \code{obj} is a path/filename.
  The input data file if specified must be a plain text file with variable names in the first row as readable via the command \code{read.table(datafilename, header = TRUE)}.
}
\value{%
  \code{pattLrep.fit} returns an object of class \code{pattMod}.
  The function \code{\link{print}} (i.e., \code{\link{print.pattMod}}) can be used to print the results and the function \code{\link{patt.worth}} to produce a matrix of worth parameters.
  An object of class \code{pattMod} is a list containing  the following components:

  \item{}{main results of the fit like estimates (\code{coefficients}), log likelihood of the model (\code{ll}), log likelihood of the saturated model (\code{fl}), and the call}
  \item{result}{a list of results from the fitting routine (see Value of \code{\link{nlm}}.}
  \item{envList}{a list with further fit details like subject covariates design structure \code{covdesmat}, paired comparison response pattern matrix \code{Y}, etc.}
  \item{partsList}{%
    a list of the basic data structures for each subgroup defined by crossing all covariate levels and different missing value patterns.
    Each element of \code{partsList} is again a list containing counts, missing value pattern, the CL matrix represented as the vector \code{s}, and the specification of the covariates.
    Use \code{\link{str}} to inspect the elements and see example below.
  }
}
\note{%
  The number of timepoints can also be regarded as different response dimensions.
}
\author{%
  Reinhold Hatzinger
}
\section{Warning}{%
  The size of the table to be analysed increases dramatically with the number of items and time points.
  For ratings (Likert items) the number of paired comparison response categories is always three.
  For each time point the number of rows of the table to set up the design matrix is initially \eqn{(2 * number of categories - 1) ^ (number of items)}.
  After reducing to three categories the number of patterns are 13, 75, 541 for 3 to 5 items, respectively.
  Generally, the number of rows in the design matrix is \eqn{(number of patterns) ^ (number of timepoints)}.
  The number of covariate levels and the number of missing value patterns have effects only on the run time.
  A (reasonable) maximum number of items for two time points is 5, for three timepoints 4, and for four to five timepoints 3.
}
\seealso{%
  \code{\link{pattL.fit}},
  \code{\link{patt.design}},
  \code{\link{pattPC.fit}},
  \code{\link{pattR.fit}},
  \code{\link{pattRrep.fit}}
}
\examples{
# simulated data: 3 items, 2 timepoints
dat <- as.data.frame(matrix(sample(1:5, 300, replace = TRUE), ncol = 6))
res <- pattLrep.fit(dat, nitems = 3, tpoints = 2, iaT = TRUE)
res
patt.worth(res, obj.names = LETTERS[1:3])
}
\keyword{models}
\keyword{multivariate}
