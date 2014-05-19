\encoding{UTF-8}
\name{llbtPC.fit}
\alias{llbtPC.fit}
\title{Function to fit a loglinear Bradley-Terry model for paired comparisons}
\description{Function to fit a loglinear Bradley-Terry for paired comparisons
       allowing subject covariates and undecided response categories.}
\usage{
llbtPC.fit(obj, nitems, formel = ~1, elim = ~1, resptype = "paircomp",
obj.names = NULL, undec = TRUE)
}
\arguments{
  \item{obj}{either a dataframe or the path/name of the datafile to be read.
        }
  \item{nitems}{the number of compared objects, not the number of comparisons}
  \item{formel}{the formula for subject covariates to fit different
        preference scales for the objects (see below).}
  \item{elim}{the formula for the subject covariates that specify the table
        to be analysed. If omitted and \code{formel} is not \code{~1} then
        \code{elim} will be set to the highest interaction between all terms
        contained in \code{formel}. If \code{elim} is specified, the terms
        must be separated by the \code{*} operator.}
  \item{resptype}{is \code{"paircomp"} by default and is reserved for future usage.
        Any other specification will not change the behaviour of \code{llbtPC.fit}
        }
  \item{obj.names}{character vector with names for objects.
        }
  \item{undec}{for paired comparisons with a undecided/neutral category,
        a common parameter will be estimated if \code{undec = TRUE}.}
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
        allows for the successive fitting of nested models to enable the use of
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
 \code{llbtPC.fit} returns an object of class \code{llbtMod}. This object
  is basically a \code{gnm} object with an additional element \code{envList}.
  This is a list with further details like the subject covariates
  design structure \code{covdesmat}, the model specification (\code{formel}
and \code{elim}), the object names (\code{obj.names}), the number of
items (\code{nobj}) and comparisons (\code{ncomp}), etc.

The function \code{\link{llbt.worth}} can be used to
       produce a matrix of estimated worth parameters.
}
%\references{}
\author{Reinhold Hatzinger}
\note{The function \code{llbtPC.fit} is a wrapper function for \code{gnm} and
was designed to facilitate fitting of LLBTs with subject covariates and
undecided categories. More specialised setups (e.g., object-specific covariates)
can be obtained using
\code{\link{llbt.design}} and then calling \code{gnm} (or \code{\link{glm}})
directly (see Examples for \code{\link{llbt.design}}).


}
\seealso{\code{\link{llbt.design}}, \code{\link{pattL.fit}}}
\examples{
# cems universities example
res0 <- llbtPC.fit(cemspc, nitems = 6, formel = ~1,   elim = ~ENG, undec = TRUE)
res1 <- llbtPC.fit(cemspc, nitems = 6, formel = ~ENG, elim = ~ENG, undec = TRUE)

anova(res1, res0)
llbt.worth(res1)
}
\keyword{models}
\keyword{multivariate}
