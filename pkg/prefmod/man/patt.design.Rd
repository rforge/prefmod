\name{patt.design}
\alias{patt.design}
\title{
    Paired Comparison Patterns - Design Matrix Generation
}
\description{
The function \code{patt.design} converts (i) real paired comparison
responses, or (ii) a set of ratings (or Likert type responses measured on a
commmon scale), or (iii) full rankings into paired comparison patterns,
returning a new data frame containing the design matrix for a loglinear
paired comparison model.  Additionally, the frequencies of these
patterns are computed and are stored in the first column of the data
frame.
}
\usage{
patt.design(obj, nitems = NULL, objnames = "", objcovs = NULL,
            cat.scovs = NULL, num.scovs=NULL, resptype = "paircomp",
            reverse = FALSE, ia = FALSE, casewise = FALSE, ...)
}
\arguments{
  \item{obj}{either a data frame or a data file name.}
  \item{nitems}{number of items (objects). \code{nitems} is not the
           number of comparisons!
       }
  \item{objnames}{an optional character vector with names for the objects
           These names are the columns names in the ouput data frame.
           If \code{objnames} is not specified \code{o1},\code{o2}, etc.
           will be used.
       }
  \item{objcovs}{an optional data frame with object specific covariates. The rows correspond
           to the objects, the columns define the covariates. The column names of this data frame
           are later used to fit the covariates. Factors are not allowed. In that case dummy
           variables have to be set up manually (favourably using \code{\link{model.matrix}}).
       }
  \item{cat.scovs}{a character vector with the names of the categorical subject covariates
           in the data file to be included into the design matrix.
           (example: \code{cat.scovs = c("SEX", "WORK")}). If all covariates in the data
           are categorical and should be included,  the specification can be abbreviated to
           \code{cat.scovs = "ALL"}. In that case, \code{num.scovs} must not be specified.
           For no categorical covariates: \code{cat.scovs = ""}, the default.
       }
  \item{num.scovs}{analogous to  \code{cat.scovs} for numerical (continuous) subject covariates. If
           any numerical covariates are specified, \code{casewise} is set to \code{TRUE}.
       }
  \item{resptype}{one of \code{"paircomp"}, \code{"rating"}, or \code{"ranking"}.
       }
  \item{reverse}{If the responses are such that
           low values correspond to high preference (or agreement or rank) and high values
           to low preference (or agreement or ranks) (e.g., (1) \emph{I strongly agree} \ldots
           (5)  \emph{I strongly disagree}) then \code{blnRevert}
           should be specified to be \code{FALSE}, the default. Otherwise
           set \code{reverse = TRUE}. The only exception is paired comparison responses
           that are coded -1/1, -1/0/1, -2/-1/0/1/2, etc. Then negative numbers are treated as not preferred.
           (See Input Data below)
       }
  \item{ia}{generates covariates for interactions between comparisons if
           \code{ia = TRUE}.
       }
  \item{casewise}{If \code{casewise = TRUE} a separate design structure is set up
       for each subject in the data. This is required when fitting continuous
       subject covariates. However, the design can become very large in the case of many subjects
       and/or comparisons. See Details below.
       }
  \item{...}{deprecated options to allow for backwards compatibility (see Deprecated below)
  }

}
\details{


   The function \code{patt.design} allows for different scenarios mainly concerning
   \itemize{
       \item \strong{responses.} Currently, three types of responses can be specified.
             \itemize{
                 \item \strong{\emph{paired comparison data.}} Responses can be either simply
                    \emph{preferred} -- \emph{not preferred} or ordinal
                    (\emph{strongly preferred} -- \dots -- \emph{not at all preferred}).
                    In both cases an undecided category may or may not occur. If there are more
                    than three categories a they are reduced to two or three response categories.
                    The set of paired
                    comparison responses represents a response pattern.
                 \item \strong{\emph{ratings/Likert type responses.}}
                    The responses to Likert type items are transformed
                    to paired comparison responses by calculating the difference
                    between each pair of the Likert items. This leads to an ordinal (adjacent
                    categories) paired comparison model with 2\eqn{k}-1 response categories where
                    \eqn{k} is the number
                    of the (original) Likert categories. Again, the transformed ratings are
                    reduced to three response categories
                    (\emph{preferred} -- \emph{undecided} -- \emph{not preferred}).
                 \item \strong{\emph{rankings.}} Currently only full rankings are allowed, i.e.,
                    a (consecutive) integer must uniquely be assigned to each object in a list
                    according to the (subjective) ordering. Ties are not allowed.
                    As for ratings,
                    the rankings are transformed
                    to paired comparison responses by calculating the difference
                    between each pair of the ranks. Again a category reduction (as described above) is
                    automatically performed.
             }
       \item \strong{comparison covariates.} The design matrix for the basic model has columns
             for the items (objects) and (depending on the type of responses)
             for undecided comparisons. For ratings (Likert type) undecided comparisons occur
             if any subject has responded to two items in the same category. For paired comparisons
             it depends on the design. For rankings there are no undecided categories.
             If undecided categories occur
             there is one dummy variable for each comparison.
             Additionally, covariates for two way interaction between comparisons
             (i.e., for effects resulting from the dependence between two
             comparisons that have one item in common) can be obtained by setting
             \code{ia = TRUE}.
       \item \strong{object specific covariates.} For modelling certain characteristics
             of objects a reparameterisation can be included in the design. This is sometimes
             called conjoint analysis. The object specific covariates can be continuous or
             dummy variables. For the specification see Argument \code{objcovs} above.
       \item \strong{subject covariates.} For modelling different preference
             scales for the items according to characteristics of the respondents
             categorical subject covariates can be included in the design.
             The corresponding variables are defined as numerical vectors where
             the levels are specified with consecutive integers starting with 1.
             This format must be used in the input data file and is also used in
             all outputs.
   }
}

\section{Input Data}{
        Responses have to be coded as consecutive integers (e.g.,
        (0,1), or (1,2,3,\dots), where the smallest value corresponds to
        (highest) preference for the first object in a comparison. 
	
        For (ordinal) paired comparison data (\code{resptype = "paircomp"}) the codings
        (1,-1), (2,1,-1,-2), (1,0,-1), (2,1,0,-1,-2) etc. can also be used.
        Then negative numbers correspond to not preferred, 0 to undecided.
	Missing responses are not allowed (use functions
        \code{\link{pattPC.fit}}, \code{\link{pattL.fit}}, or \code{\link{pattR.fit}} instead). 

        Input data (via the first argument \code{obj} in the function call)
        is specified either through a dataframe or
        a datafile in which case \code{obj} is a path/filename.  The input
        data file if specified must be a plain text file with variable names in
        the first row as readable via the command \code{read.table(datafilename,
        header = TRUE)}.

        The leftmost columns must be the
        responses to the paired comparisons, ratings (Likert items), or rankings.
        For paired comparisons the mandatory order is of comparisons is
        (12) (13) (23) (14) (24) (34) (15) (25) etc. For rankings, the lowest value
        means highest rank according to the underlying scale.
        Each column in the data file corresponds to one of the ranked objects. For example,
        if we have 3 objects denoted by \code{A},\code{B},and \code{C}, with
        corresponding columns in the data matrix, the response pattern \code{(3,1,2)}
        represents: object \code{B} ranked highest, \code{C} ranked second, and
        \code{A} ranked lowest. For ratings. again the lowest value means highest
        \sQuote{endorsement} (agreement) according to the underlying scale. All items
        are assumed to have the same number of response cateory

        The columns for responses are optionally followed  by columns for 
	subject covariates. If categorical, they have to be specified such that the categories are represented
        by consecutive integers starting with 1. Missing values are not allowed 
        and treated such that
        rows with NAs are removed from the resulting design structure and a message is printed.
        For an example see \code{\link{xmpl}}.

        (Besides supplying data via a dataframe or a datafile name, \code{obj} can also
        be specified as a control list with the same elements as the arguments in the function call.
        The data must then be specified as a path/filename using the element \code{datafile = "filename"}.
        The control list feature is deprecated. An example is given below.)

}

\value{
       The output is a dataframe. Each row represents a unique response pattern.
       If subject covariates are specified, each row instead represents a particular
       combination of a unique covariate combination with a response pattern. All
       possible combinations are generated.\cr
       The first column contains the counts for the paired
       comparison response patterns and is labelled with \code{Y}. The next columns
       are the covariates for the items and the undecided category effects (one for
       each comparison). These are labelled as \code{u12},\code{u13}, etc., where
       \code{12} denotes the
       comparison between items \code{1} and \code{2}. Optionally,
       covariates for dependencies between comparisons follow. The columns are labelled
       \code{Ia.bc} denoting the interaction of the comparisons between items \code{(a,b)}
       and \code{(a,c)} where the common item is \code{a}. If subject covariates are
       present they are in the rightmost columns and defined to be factors.\cr
%%       The response patterns themselves are not given in the output dataframe.
%%       However,
%%       the rows of the dataframe are in standard order, with the patterns varying within
%%       covariate combinations, and the patterns varying with the leftmost comparison
%%       (1,2), (1,3), (2,3), (1,4) etc varying fastest.
}
\section{Deprecated }{The following options are for backwards compatibility and should
       no longer be used..
\describe{
  \item{blnCasewise}{same as \code{casewise}.}
  \item{cov.sel}{same as \code{cat.scovs}.}
}
Options for requesting GLIM commands and data structures are no longer supported.
Specifying the input to \code{llbt.design} via a control list is also deprecated.
If you want to use these features you have to install prefmod <= 0.8-22.
}
\references{Dittrich, R., Francis, B.J., Hatzinger R., Katzenbeisser, W. (2007),
A Paired Comparison Approach for the Analysis of Sets of Likert Scale Responses.
\emph{Statistical Modelling}, Vol. 7, No. 1, 3-28. }
\author{Reinhold Hatzinger}
\seealso{\code{\link{llbt.design}},
\code{\link{pattPC.fit}}, \code{\link{pattL.fit}}, \code{\link{pattR.fit}}
}
\examples{

## mini example with three Likert items
## and two subject covariates

data(xmpl)   # example data in package
dsgnmat <- patt.design(xmpl, nitems = 3, resptype="rating",
    blnIntcovs = TRUE, cov.sel="ALL")
print(head(dsgnmat))


## ILLUSTRATING THE ISSP2000 EXAMPLE
##   simplified version of the analysis as
##   given in Dittrich et.al.(2007)

data(issp2000)

design <- patt.design(issp2000, nitems=6, resptype="rating",
      cov.sel=c("SEX","EDU"))


#  - fit null multinomial model (basic model for items without
#    subject covariates) through Poisson distribution.
#  - SEX:EDU parameters are nuisance parameters
#  - the last item (GENE) becomes a reference item
#    in the model and is aliased; all other items
#    are compared to this last item

# item parameters with undecided effects and no covariate effects.

summary(glm(y~SEX*EDU + CAR+IND+FARM+WATER+TEMP+GENE
                + u12+u13+u23+u14+u24+u34+u15+u25+u35+u45+u16+u26+u36+u46+u56,
             family=poisson, data=design))

# now add main effect of SEX on items

summary(glm(y~SEX:EDU + CAR+IND+FARM+WATER+TEMP+GENE
                + (CAR+IND+FARM+WATER+TEMP+GENE):SEX
                + u12+u13+u23+u14+u24+u34+u15+u25+u35+u45+u16+u26+u36+u46+u56,
             family=poisson, data=design))
}
\keyword{models}
\keyword{regression}
