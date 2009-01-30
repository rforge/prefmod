\name{llbt.design}
\alias{llbt.design}
\title{Loglinear Bradley-Terry Model (LLBT) - Design Matrix Generation}
\description{
The function \code{llbt.design}
returns a data frame containing the design matrix for a loglinear
paired comparison model.  Additionally, the frequencies of the
pairwise comparisons are computed and are stored in the first column of the data
frame.  Optionally, the function provides all necessary structures
(commands, data/design files) to fit the loglinear paired comparisons
pattern model in GLIM, which is often more efficient at fitting large
loglinear models of this type.
}
\usage{
llbt.design(obj, nitems = NULL, objnames = "",
    blnCasewise = FALSE, cov.sel = "",
    blnGLIMcmds = FALSE, glimCmdFile = "", outFile = "")
}
\arguments{
  \item{obj}{either a data frame, a data file name, or a
           control object.
       }
  \item{nitems}{number of items (objects).
       }
  \item{objnames}{an optional character vector with names for the objects
           These names are the columns names in the ouput data frame.
           If \code{objnames} is not specified \code{o1},\code{o2}, etc.
           will be used.
       }
  \item{blnCasewise}{If \code{blnCasewise = TRUE} a separate design structure is set up
       for each subject in the data. This is required when fitting continuous
       subject covariates. However, the design can become very large in the case of many subjects
       and/or comparisons. See Details below.}
  \item{cov.sel}{a character vector with the names of the subject covariates
           in the data file to be included into the design matrix.
           (example: \code{cov.sel = c("SEX", "AGE")}). If all covariates are to be
           included the specification can be abbreviated to \code{cov.sel = "ALL"}.
           For no covariates specify \code{cov.sel = ""}.
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
}
\details{
   The function \code{llbt.design} allows for different scenarios mainly concerning
   \itemize{
      \item \strong{\emph{paired comparison data.}} Responses can be either simply
         \emph{preferred} -- \emph{not preferred} or ordinal
         (\emph{strongly preferred} -- \dots -- \\emph{not at all preferred}).
         In both cases an undecided category may or may not occur. If there are more
         than three categories a they are reduced to two or three response categories.
      \item \strong{item covariates.} The design matrix for the basic model has columns
         for the items (objects) and for each response category
      \item \strong{subject covariates.} For modelling different preference
         scales for the items according to characteristics of the respondents
         categorical subject covariates can be included in the design.
         The corresponding variables are defined as numerical vectors where
         the levels are specified with consecutive integers starting with 1.
         This format must be used in the input data file and is also used in
         all outputs. Prior to fit a model, the categorical subject covariates have to
         be declared as \code{\link{factor}(s)}. See example below.
       \item \strong{GLIM output}. If the user specifies \code{blnGLIMcmds = TRUE}
             two files are generated one of which contains all GLIM commands to fit
             a basic loglinear paired comparisons pattern model. The other contains
             the design matrix optionally including subject covariates.
   }
}
\section{Input Data}{
        Responses have to be coded as consecutive integers (e.g.,
        (0,1), or (1,2,3,\dots), where the lowest value corresponds to
        (highest) preference for the first object in a comparison. For paired comparison without undecided category
        (-1,1), or (1,0,-1) for paired comparison with an undecided category, can also be used
        ('-1' is the not preferred category). Missing responses (for paired comparisons but not for subject covariates)
        are allowed
        under a missing at random assumption and specified via \code{NA}.

        Input data (via the first argument \code{obj} in the function call)
        is specified either through a dataframe or
        a datafile in which case \code{obj} is a path/filename.  The input
        data file if specified must be a plain text file with variable names in
        the first row as readable via the command \code{read.table(datafilename,
        header = TRUE)}.

        The leftmost columns must be the
        responses to the paired comparisons (where the
        mandatory order of comparisons is (12) (13) (23) (14) (24) (34) (15)
        (25) etc.),
        optionally followed  by columns for
        subject covariates. If categorical, these have to be specified such that the categories are represented
        by consecutive integers starting with 1. Missing values for subject covariates are not allowed
        and treated such that
        rows with NAs are removed from the resulting design structure and a message is printed.


        For an example see \code{\link{xmpl}} or the file \code{xmpl.dat} in the package's
        \code{data/} directory.

        (Besides supplying data via a dataframe or a datafile name, \code{obj} can also
        be specified as a control list with the same elements as the arguments in the function call.
        The data must then be specified as a path/filename using the element \code{datafile = "filename"}.
        The control list feature is deprecated. An example with slightly different specification
        is given in \code{\link{patt.design}}.)

}
\value{
       The output is a dataframe. Each row represents a decision in a certain comparison.
       Dependent on the number of response categories, comparisons are made up of two or
       three rows in the design matrix.
       If subject covariates are specified, the design matrix is duplicated as many times as there
       are combinations of the levels of each categorical covariate or, if
       \code{blnCasewise = TRUE}, as there are subjects in the data. Each individual
       design matrix consists of rows for all comparisons.

       The first column contains the counts for the paired
       comparison response patterns and is labelled with \code{y}. The next columns
       are the covariates for the categories (labelled as \code{g0},\code{g1}, etc.)
       and for the items. If subject covariates are
       present they are in the rightmost columns.

       Alternatively, the function \code{llbt.design} does not produce visible output in \R if GLIM output
       is requested via \code{blnGLIMcmds = TRUE}. The output is then written to
       the corresponding files.
}
\references{R. Dittrich, R. Hatzinger and W. Katzenbeisser.
Modelling the effect of subject-specific covariates in paired
comparison studies with an application to university rankings.
Applied Statististics (1998), 47, Part 4,  pp. 511-525
}
\author{Reinhold Hatzinger}
\seealso{\code{\link{patt.design}},
\code{\link{pattPC.fit}}, \code{\link{pattL.fit}}, \code{\link{pattR.fit}}
}
\examples{
## cems universities example
data(cemspc)
des<-llbt.design(cemspc,nitems=6,cov.sel="ENG")

eng<-factor(des$ENG)
res<-glm(y~o1+o2+o3+o4+o5+o6+eng:(o1+o2+o3+o4+o5+o6)+mu*eng,
      family=poisson, data=des)
summary(res)

## calculating and plotting worth parameters
lambda<-res$coefficients[c(2:6,23:27)]
lambda[6:10]<-lambda[1:5]+lambda[6:10]
lambmat<-matrix(lambda,ncol=2)
lambmat<-rbind(lambmat,c(0,0))
worthmat<-apply(lambmat,2,
      function(x) exp(2*x)/sum(exp(2*x)))
colnames(worthmat)<-c("English:good","English:poor")
rownames(worthmat)<-c("London", "Paris", "Milano",
      "St.Gallen", "Barcelona","Stockholm")
plotworth(worthmat)
}
\keyword{models}
\keyword{regression}
