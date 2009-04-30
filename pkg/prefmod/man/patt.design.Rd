\name{patt.design}
\alias{patt.design}
\title{
    Paired Comparison Patterns - Design Matrix Generation
}
\description{
The function \code{patt.design} converts (i) real paired comparison
responses, or (ii) a set of ratings (Likert type responses measured on a
commmon scale), or (iii) full rankings into paired comparison patterns,
returning a new data frame containing the design matrix for a loglinear
paired comparison model.  Additionally, the frequencies of these
patterns are computed and are stored in the first column of the data
frame.  Optionally, the function provides all necessary structures
(commands, data/design files) to fit the loglinear paired comparisons
pattern model in GLIM, which is often more efficient at fitting large
loglinear models of this type.
}
\usage{
patt.design(obj, nitems = NULL, objnames = "", resptype = "paircomp",
     blnRevert = FALSE, cov.sel = "", blnIntcovs = FALSE,
     blnGLIMcmds = FALSE, glimCmdFile = "", outFile = "", intFile = "")
}
\arguments{
  \item{obj}{either a data frame, a data file name, or a
           control object.
       }
  \item{nitems}{number of items (objects). \code{nitems} is not the
           number of comparisons!
       }
  \item{objnames}{an optional character vector with names for the objects
           These names are the columns names in the ouput data frame.
           If \code{objnames} is not specified consecutive letters
           starting with "a" will be used.
       }
  \item{resptype}{one of \code{"paircomp"}, \code{"rating"}, or \code{"ranking"}.
       }
  \item{blnRevert}{If the responses are such that
           low values correspond to high preference or agreement and high values
           to low preference or agreement (e.g., (1) \emph{I strongly agree} \ldots
           (5)  \emph{I strongly disagree}) then \code{blnRevert}
           should be specified to be \code{FALSE}. Otherwise
           set \code{blnRevert = TRUE}.
       }
  \item{cov.sel}{a character vector with the names of the subject covariates
           in the data file to be included into the design matrix.
           (example: \code{cov.sel = c("SEX", "AGE")}). If all covariates are to be
           included the specification can be abbreviated to \code{cov.sel = "ALL"}.
           For no covariates specify \code{cov.sel = ""}.
       }
  \item{blnIntcovs}{generates covariates for interactions between comparisons if
           \code{blnIntcovs = TRUE}.
       }
  \item{blnGLIMcmds}{\code{TRUE}, if GLIM output is wanted. If \code{blnGLIMcmds = FALSE}
           the following items can be set to any value (such as a null text string) and are ignored.
           Please note that if \code{blnGLIMcmds} is set to be \code{TRUE} there is no output in \R
           but instead goes to the the following files.
       }
  \item{glimCmdFile}{name of the output file which will contain all necessary commands to
          fit a basic model (defining all structures and reading the
          necessary data).
       }
  \item{outFile}{name of the data/design file to be read into GLIM. It consists
           of the response frequencies and the covariates for the objects,
           the undecided comparison responses and the subject effects.
       }
  \item{intFile}{name of the design file for the interaction effects.
          This file is only generated if \code{blnIntcovs = TRUE}.
       }
}
\details{


   The function \code{patt.design} allows for different scenarios mainly concerning
   \itemize{
       \item \strong{responses.} Currently, three types of responses can be specified.
             \itemize{
                 \item \strong{\emph{paired comparison data.}} Responses can be either simply
                    \emph{preferred} -- \emph{not preferred} or ordinal
                    (\emph{strongly preferred} -- \dots -- \\emph{not at all preferred}).
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
       \item \strong{item covariates.} The design matrix for the basic model has columns
             for the items (objects) and (depending on the type of responses)
             for undecided comparisons. For ratings (Likert type) undecided comparisons occur
             if any subject has responded to two items in the same category. For paired comparisons
             it depends on the design. For rankings there are no undecided categories.
             If undecided categories occur
             there is one dummy variable for each comparison.
             Additionally, covariates for two way interaction between comparisons
             (i.e., for effects resulting from the dependency between two
             comparisons that have one item in common) can be obtained by setting
             \code{blnIntcovs = TRUE}.
       \item \strong{subject covariates.} For modelling different preference
             scales for the items according to characteristics of the respondents
             categorical subject covariates can be included in the design.
             The corresponding variables are defined as numerical vectors where
             the levels are specified with consecutive integers starting with 1.
             This format must be used in the input data file and is also used in
             all outputs.
       \item \strong{GLIM output}. If the user specifies \code{blnGLIMcmds = TRUE}
             two files are generated one of which contains all GLIM commands to fit
             a basic loglinear paired comparisons pattern model. The other contains
             the design matrix optionally including subject covariates. If
             dependency covariates are requested they are written to a third file.
             (Please note that the corresponding part of design matrix is transposed
             in the interactions output file to allow GLIM for using the \code{$array}
             facility in case of a large number of parameters to be estimated)
   }
}

\section{Input Data}{
        Responses have to be coded as consecutive integers (e.g.,
        (0,1), or (1,2,3,\dots). For paired comparison without undecided category
        (-1,1), or (1,0,-1) for paired comparison with an undecided category, can also be used
        ('-1' is the not preferred category).

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

        The columns for responses are optionally followed  by columns for categorical
        subject covariates. These have to be specified such that the categories are represented
        by consecutive integers starting with 1. Missing values are not allowed (use functions
        \code{\link{pattPC.fit}}, \code{\link{pattL.fit}}, or \code{\link{pattR.fit}} instead),
        and treated such that
        rows with NAs are removed from the resulting design structure and a message is printed.
        For an example see \code{\link{xmpl}} or the file \code{xmpl.dat} in the package's
        \code{data/} directory.

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
       Alternatively, the function \code{patt.design} does not produce output in \R if GLIM output
       is requested via \code{blnGLIMcmds = TRUE}. The output is then written to
       the corresponding files (see \bold{The Control List} below).
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

\dontrun{
## example with control list
## not run because input data file does not exist


# defining the ctrl list
# would be typically read from file using source()
testex1<-list(
   resptype     = "paircomp",
   datafile     = "test/test.dat",
   nitems       = 5,
   blnRevert    = FALSE,
   blnReducecat = TRUE,
   blnIntcovs   = FALSE,
   cov.sel      = c("SEX","URB"),

   blnGLIMcmds  = TRUE,
   glimCmdFile  = "test/test.gli",
   outFile      = "test/test.design",
   intFile      = ""        # since blnIntcovs = FALSE
)
patt.design(testex1)
}


}
\keyword{models}
\keyword{regression}
