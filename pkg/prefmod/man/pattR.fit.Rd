\name{pattR.fit}
\alias{pattR.fit}
\title{Function to fit a pattern model for (partial) rankings}
\description{Function to fit a pattern model for (partial) rankings (transformed
      to paired comparisons)
       allowing for missing values using a CL approach.}
\usage{
pattR.fit(obj, nitems, formel = ~1, elim = ~1, resptype = "ranking",
    ia = FALSE, NItest = FALSE, pr.it = FALSE)
}
\arguments{
  \item{obj}{either a dataframe or the path/name of the datafile to be read.
       }
  \item{nitems}{the number of items}
  \item{formel}{the formula for subject covariates to fit different
        preference scales for the objects (see below).}
  \item{elim}{the formula for the subject covariates that specify the table
        to be analysed. If ommitted and \code{formel} is not \code{ ~1} then
        \code{elim} will be set to the highest interaction between all terms
        contained in \code{formel}. If \code{elim} is specified, the terms
        must be separated by the \code{*} operator.}
  \item{resptype}{is \code{"ranking"} by default and is reserved for future usage.
        Any other specification will not change the behaviour of \code{pattR.fit}
        }
  \item{ia}{interaction parameters between comparisons that have one object
        in common if \code{ia = TRUE}.
        }
  \item{NItest}{separate estimation of object parameters for complete and
        incomplete patterns if \code{NItest = TRUE}. Currently,
        \code{NItest} is set to \code{FALSE} if subject covariates are specified.
        }
  \item{pr.it}{a dot is printed at each iteration cycle if set to \code{TRUE}
        }
}
\details{Models including categorical subject covariates can be fitted using the
        \code{formel} and \code{elim} arguments. \code{formel} specifies the
        actual model to be fitted. For instance, if specified as
        \code{formel=~SEX} different preference scale for the objects will be
        estimated for males and females. For two or more covariates,
        the operators \code{+} or \code{*} can be used to model main or interaction
        effects, respectively. The operator \code{:} is not allowed. See also
        \code{\link{formula}}.

        The specification for \code{elim} follows the same rules as for
        \code{formel}.  However, \code{elim} specifies the basic contingency
        table to be set up but does not specify any covariates to be fitted.
        This is done using \code{formel}.
        If, e.g., \code{elim=~SEX} but \code{formel=~1},
        then the table is set up as if \code{SEX} would be fitted but only one global
        preference scale is computed.  This feature
        allows for the succesive fitting of nested models to enable the use of
        deviance differences for model selection (see example below).
}
\section{Input Data}{
       The responses have to be coded as consecutive integers starting with 1.
       The value of 1 means highest rank according to the underlying scale.
       Each column in the data file corresponds to one of the ranked objects. For example,
       if we have 3 objects denoted by \code{A},\code{B},and \code{C}, with
       corresponding columns in the data matrix, the response pattern \code{(3,1,2)}
       represents: object \code{B} ranked highest, \code{C} ranked second, and
       \code{A} ranked lowest. Missing values are coded as \code{NA},
       ties are not allowed (in that case use \code{\link{pattL.fit}}.
       Rows with less than 2 ranked objects are removed from the fit
       and a message is printed.

       Optional subject covariates have to be specified
       such that the categories are represented by consecutive integers
       starting with 1. Rows with missing values for subject covariates are
       removed from the data and a message is printed.  The leftmost columns in
       the data must be the rankings, optionally followed by columns for
       categorical subject covariates.

       The data specified via \code{obj} are supplied using either a data frame
       or a datafile in which case \code{obj} is a path/filename.  The input
       data file if specified must be a plain text file with variable names in
       the first row as readable via the command \code{read.table(datafilename,
       header = TRUE)}.

       For an example without covariates and no mising values
       see \code{\link{salad}}.
       }
\value{
  \code{pattR.fit} returns an object of class \code{pattMod}. The function
       \code{\link{print}} (i.e., \code{\link{print.pattMod}}) can be used to
       print the results and the function \code{\link{patt.worth}} to
       produce a matrix of worth parameters. An object of class \code{pattMod}
       is a list containing  the following components:
     \item{}{main results of the fit like estimates (\code{coefficients}),
     log likelihood of the model (\code{ll}), log likelihood of the saturated model
     (\code{fl}), and the call}
     \item{result}{a list of results from the fitting routine (see Value of
          \code{\link{nlm}}).
          }
     \item{envList}{a list with further fit details like subject covariates
          design structure \code{covdesmat}, paired comparison reponse pattern matrix
          \code{Y}, etc.}
     \item{partsList}{a list of the basic data structures for each subgroup
          defined by crossing all covariate levels and different missing value
          patterns. Each element of \code{partsList} is again a list containing
          counts, missing value pattern, the CL matrix represented as a vector, and
          the specification of the covariates. Use \code{\link{str}} to inspect
          the elements and see example below.
          }
      }
%\references{}
\author{Reinhold Hatzinger}
%\note{}
\section{Warning}{The size of the table to be analysed increases dramatically
     with the number of items . For rankings the number of
     paired comparison response categories is always two. The number of
     rows of the table used to set up the design matrix is \code{factorial(number of items)}.
     For instance, for nine objects this is 362880.
     A reasonale maximum number of items is 8.

     The option \code{NItest = TRUE} has to be used with care. The meaning of missing
     responses is not obvious with partial rankings. Are the corresponding values
     really missing or just not chosen.
     }
\seealso{\code{\link{patt.design}}, \code{\link{pattL.fit}}, \code{\link{pattPC.fit}}
     }
\examples{
## fit of Critchlov & Fligner (1991) Salad Dressings Data

data(salad)
pattR.fit(salad, nitems=4)

# alternatively use glm() with patt.design()

sal<-patt.design(salad,nitems=4,resptype="ranking")
glm(y~A+B+C+D,family=poisson,data=sal)
}
%
%   # Example for object covariates
%   # fit object covariates:
%   # salads A - D have varying concentrations of acetic and gluconic acid.
%   # The four pairs of concentrations are
%   # A = (.5, 0), B = (.5, 10.0), C = (1.0, 0), and D = (0, 10.0),
%
%   conc<-matrix(c(.5,0, .5,10,  1,0, 0,10), nc=2, byrow=TRUE)
%   X<-as.matrix(sal[,2:5])%*%conc
%   glm(y~X,family=poisson,data=sal)
%
\keyword{models}
\keyword{multivariate}
