\name{pattPC.fit}
\alias{pattPC.fit}
\title{Function to fit a pattern model for paired comparisons}
\description{Function to fit a pattern model for paired comparisons
       allowing for missing values using a CL approach.}
\usage{
pattPC.fit(obj, nitems, formel = ~1, elim = ~1, resptype = "paircomp",
    obj.names = NULL, undec = FALSE, ia = FALSE, NItest = FALSE,
    NI = FALSE, MIScommon = FALSE, MISmodel="obj", MISalpha = NULL,
    MISbeta = NULL, pr.it = FALSE)
}
\arguments{
  \item{obj}{either a dataframe or the path/name of the datafile to be read.
        }
  \item{nitems}{the number of compared objects, not the number of comparisons}
  \item{formel}{the formula for subject covariates to fit different
        preference scales for the objects (see below).}
  \item{elim}{the formula for the subject covariates that specify the table
        to be analysed. If ommitted and \code{formel} is not \code{ ~1} then
        \code{elim} will be set to the highest interaction between all terms
        contained in \code{formel}. If \code{elim} is specified, the terms
        must be separated by the \code{*} operator.}
  \item{resptype}{is \code{"paircomp"} by default and is reserved for future usage.
        Any other specification will not change the behaviour of \code{pattPC.fit}
        }
  \item{obj.names}{character vector with names for objects.
        }
  \item{undec}{for paired comparisons with a undecided/neutral category,
        a common parameter will be estimated if \code{undec = TRUE}.}
  \item{ia}{interaction parameters between comparisons that have one object
        in common if \code{ia = TRUE}.
        }
  \item{NItest}{separate estimation of object parameters for complete and
        incomplete patterns if \code{NItest = TRUE}. Currently,
        \code{NItest} is set to \code{FALSE} if subject covariates are specified.
        }
  \item{NI}{if  \code{TRUE}, fits large table (crossclassification with NA patterns), for comparison
        with models including \code{MISalpha} (and \code{MISbeta}). }
  \item{MIScommon}{if \code{TRUE}, fits a common parameter for NA indicators,
        i.e.,
        \eqn{\alpha=\alpha_i=\alpha_j=\dots}{alpha = alpha_i = alpha_j = ...}.}
  \item{MISmodel}{either \code{obj} (default) or \code{comp}. Specifies the
        model for missing values. If \code{obj}, parameters for missing objects
        or if \code{comp}, parameters for missing comparisons are estimated. In both
        cases \code{MISalpha} (and optionally \code{MISbeta}) specify which parameters
        are involved.
        }
  \item{MISalpha}{if not \code{NULL}, specification to fit parameters for
        NA indicators using a
        logical vector, where \code{TRUE} means that the NA indicator parameter
        for the corresponding object (or comparison, dependent on \code{MISmodel})
        should be estimated (see example below).
        For \code{MISmodel="obj"} each comparison is reparameterized  with
        \eqn{\alpha_i+\alpha_j}{alpha_i + alpha_j}, for \code{MISmodel="comp"}
        each comparison is reparameterized  with \eqn{\alpha_{ij}}{\alpha_{ij}}.
        }
  \item{MISbeta}{if not \code{NULL}, fits parameters for MNAR model, i.e., interactions
        between outcome model parameters and NA indicator parameters.
        The reparameterization is analogous to \code{MISalpha}. Usually,
        the specification for \code{MISbeta} is the same as for \code{MISalpha}, but any subset is
        reasonable. If \code{MISalpha = NULL} but \code{MISbeta} is not, then
        \code{MISalpha} is set to  \code{MISbeta}.
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
       The responses have to be coded as 0/1 for paired comparisons without
       undecided category (0 means first object in a comparison preferred) or
       0/1/2 for paired comparisons with an undecided category (where 1 is the
       undecided category).  Optional subject covariates have to be specified
       such that the categories are represented by consecutive integers
       starting with 1. Rows with missing values for subject covariates are
       removed from the data and a message is printed.  The leftmost columns in
       the data must be the responses to the paired comparisons (where the
       mandatory order of comparisons is (12) (13) (23) (14) (24) (34) (15)
       (25) etc.), optionally followed by columns for
       categorical subject covariates.

       The data specified via \code{obj} are supplied using either a data frame
       or a datafile in which case \code{obj} is a path/filename.  The input
       data file if specified must be a plain text file with variable names in
       the first row as readable via the command \code{read.table(datafilename,
       header = TRUE)}.

       For an example see \code{\link{cemspc}}.
       }
\value{
  \code{pattPC.fit} returns an object of class \code{pattMod}. The function
       \code{\link{print}} (i.e., \code{\link{print.pattMod}}) can be used to
       print the results and the function \code{\link{patt.worth}} to
       produce a matrix of the estimated worth parameters. An object of class "\code{pattMod}"
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
     with the number of objects. For paired comparisons with two response
     categories the number of rows of the table is \eqn{ 2 ^ (number of comparisons)}, e.g.,
     with six objects this is 32768, for three response categories this is 14348907.
     A reasobale maximum number of objects to be analysed with pattern models is
     6 in the case of two response categories and 5 when an additional undecided/neutral
     category has been observed).
     }
\seealso{\code{\link{patt.design}}, \code{\link{checkMIS}}, \code{\link{pattL.fit}},
     \code{\link{pattR.fit}}
     }
\examples{
## fit only first three objects with undecided parameter
data(cemspc)
pattPC.fit(cemspc, nitems=3, undec=TRUE)

## check for ignorable missing
pattPC.fit(cemspc, nitems=3, undec=TRUE, NItest=TRUE)

## check if SEX has an effect
m1<-pattPC.fit(cemspc, nitems=3, formel=~1,elim=~SEX, undec=TRUE)
m2<-pattPC.fit(cemspc, nitems=3, formel=~SEX, elim=~SEX, undec=TRUE)

## calculate LR test for SEX
ll1<-m1$result$minimum
ll2<-m2$result$minimum
df1<-length(m1$result$estimate)
df2<-length(m2$result$estimate)
lr<-2*(ll1-ll2)
df<-df2-df1
cat("LR test is",lr,"on df =",df," ( p =",round(1-pchisq(lr,df),digits=4),")\n")


## generates data set with three items and some missing values
## in comparison (23), column 3, then there are no NAs for
## object 1
data(dat4)
data3<-dat4[,1:3]
idx3<-sample(1:100,10)
data3[idx3,3]<-NA
checkMIS(data3,nitems=3,verbose=TRUE)

## estimate MNAR PC pattern model for data3 without alpha1 and beta1
pattPC.fit(data3, nitems=3, MISalpha=c(FALSE,TRUE,TRUE), MISbeta=c(FALSE,TRUE,TRUE))
}
\keyword{models}
\keyword{multivariate}
