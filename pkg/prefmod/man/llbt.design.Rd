\encoding{UTF-8}
\name{llbt.design}
\alias{llbt.design}
\title{Loglinear Bradley-Terry Model (LLBT) -- Design Matrix Generation}
\description{%
  The function \code{llbt.design} returns a data frame containing the design matrix for a loglinear paired comparison model.
  Additionally, the frequencies of the pairwise comparisons are computed and are stored in the first column of the data frame.
}
\usage{%
llbt.design(data, nitems = NULL, objnames = "", objcovs = NULL,
        cat.scovs = NULL, num.scovs = NULL, casewise = FALSE, \dots)
}
\arguments{%
  \item{data}{%
    either a data frame or a data file name.
  }
  \item{nitems}{%
    number of items (objects).
  }
  \item{objnames}{%
    an optional character vector with names for the objects.
    These names are the columns names in the output data frame.
    If \code{objnames} is not specified \code{o1}, \code{o2}, etc. will be used.
  }
  \item{objcovs}{%
    an optional data frame with object specific covariates.
    The rows correspond to the objects, the columns define the covariates.
    The column names of this data frame are later used to fit the covariates.
    Factors are not allowed.
    In that case dummy variables have to be set up manually (favourably using \code{\link{model.matrix}}).
  }
  \item{cat.scovs}{%
    a character vector with the names of the categorical subject covariates in the data file to be included into the design matrix.
    (Example: \code{cat.scovs = c("SEX", "WORK")}.)
    If all covariates in the data are categorical and should be included,  the specification can be abbreviated to \code{cat.scovs = "ALL"}.
    In that case, \code{num.scovs} must not be specified.
    For no categorical covariates: \code{cat.scovs = ""}, the default.
  }
  \item{num.scovs}{%
    analogous to  \code{cat.scovs} for numerical (continuous) subject covariates.
    If any numerical covariates are specified, \code{casewise} is set to \code{TRUE}.
  }
  \item{casewise}{%
    If \code{casewise = TRUE} a separate design structure is set up for each subject in the data.
    This is required when fitting continuous subject covariates.
    However, the design can become very large in the case of many subjects and/or comparisons.
    See Details below.
  }
  \item{\dots}{%
    deprecated options to allow for backwards compatibility (see Deprecated below)
  }
}
\details{
  The function \code{llbt.design} allows for different scenarios mainly concerning
  \itemize{
    \item\strong{Paired comparison data.}
      Responses can be either simply \emph{preferred} -- \emph{not preferred} or ordinal (\emph{strongly preferred} -- \dots -- \emph{not at all preferred}).
      In both cases an undecided category may or may not occur.
      If there are more than three categories a they are reduced to two or three response categories.
    \item\strong{Item covariates.}
      The design matrix for the basic model has columns for the items (objects) and for each response category.
    \item\strong{Object specific covariates.}
      For modelling certain characteristics of objects a reparameterisation can be included in the design.
      This is sometimes called conjoint analysis.
      The object specific covariates can be continuous or dummy variables.
      For the specification see Argument \code{objcovs} above.
    \item\strong{Subject covariates.}
      For modelling different preference scales for the items according to characteristics of the respondents categorical and/or continuous subject covariates can be included in the design.\cr%
      \emph{Categorical subject covariates}: The corresponding variables are defined as numerical vectors where the levels are specified with consecutive integers starting with 1.
      This format must be used in the input data file.
      These variables are \code{\link{factor}(s)} in the output data frame.\cr%
      \emph{Continuous subject covariates}: also defined as numerical vectors in the input data frame.
      If present, the resulting design structure is automatically expanded, i.e., there are as many design blocks as there are subjects.
    \item\strong{Object specific covariates.}
      The objects (items) can be reparameterised using an object specific design matrix.
      This allows for scenarios such as conjoint analysis or for modelling some characteristics shared by the objects.
      The number of such characteristics must not exceed the number of objects minus one.
  }
}
\section{Input Data}{
  Responses have to be coded as consecutive integers (e.g., (0, 1), or (1, 2, 3, \dots), where the smallest value corresponds to (highest) preference for the first object in a comparison.
  For (ordinal) paired comparison data (\code{resptype = "paircomp"}) the codings \eqn{(1, -1)}, \eqn{(2, 1, -1, -2)}, \eqn{(1, 0, -1)}, \eqn{(2, 1, 0, -1, -2)}, etc. can also be used.
  Then negative numbers correspond to not preferred, 0 to undecided.
  Missing responses (for paired comparisons but not for subject covariates) are allowed under a missing at random assumption and specified via \code{NA}.

  Input data (via the first argument \code{obj} in the function call) is specified either through a dataframe or a datafile in which case \code{obj} is a path/filename.
  The input data file if specified must be a plain text file with variable names in the first row as readable via the command \code{read.table(datafilename, header = TRUE)}.

  The leftmost columns must be the responses to the paired comparisons (where the mandatory order of comparisons is (12) (13) (23) (14) (24) (34) (15) (25) etc.), optionally followed by columns for subject covariates.
  If categorical, these have to be specified such that the categories are represented by consecutive integers starting with 1.
  Missing values for subject covariates are not allowed and treated such that rows with NAs are removed from the resulting design structure and a message is printed.

  For an example see \code{\link{xmpl}}.
}
\value{
  The output is a dataframe of class \code{llbtdes}.
  Each row represents a decision in a certain comparison.
  Dependent on the number of response categories, comparisons are made up of two or three rows in the design matrix.
  If subject covariates are specified, the design matrix is duplicated as many times as there are combinations of the levels of each categorical covariate or, if \code{casewise = TRUE}, as there are subjects in the data.
  Each individual design matrix consists of rows for all comparisons.

  The first column contains the counts for the paired comparison response patterns and is labelled with \code{y}.
  The next columns are the covariates for the categories (labelled as \code{g0}, \code{g1}, etc.).
  In case of an odd number of categories, \code{g1} can be used to model an undecided effect.
  The subsequent columns are covariates for the items.
  If specified, subject covariates and then object specific covariates follow.

}
\section{Deprecated}{The following options are for backwards compatibility and should no longer be used.
  \describe{
    \item{blnCasewise}{same as \code{casewise}.}
    \item{cov.sel}{same as \code{cat.scovs}.}
  }
  Options for requesting GLIM commands and data structures are no longer supported.
  Specifying the input to \code{llbt.design} via a control list is also deprecated.
  If you want to use these features you have to install \pkg{prefmod} <= 0.8-22.
}
\references{%
  Dittrich, R., Hatzinger, R., & Katzenbeisser, W. (1998). Modelling the effect of subject-specific covariates in paired comparison studies with an application to university rankings.
  \emph{Journal of the Royal Statistical Society: Series C (Applied Statistics), 47}(4), 511--252.
  \doi{10.1111/1467-9876.00125}
}
\author{Reinhold Hatzinger}
\seealso{%
  \code{\link{patt.design}},
  \code{\link{llbt.worth}},
  \code{\link{llbtPC.fit}}
}
\examples{
# cems universities example
des <- llbt.design(cemspc, nitems = 6, cat.scovs = "ENG")

res0 <- gnm(y ~ o1+o2+o3+o4+o5+o6 + ENG:(o1+o2+o3+o4+o5+o6),
    eliminate = mu:ENG, data = des, family = poisson)
summary(res0)

# inclusion of g1 allows for an undecided effect
res <- gnm(y ~ o1+o2+o3+o4+o5+o6 + ENG:(o1+o2+o3+o4+o5+o6) + g1,
    eliminate = mu:ENG, data = des, family = poisson)
summary(res)

# calculating and plotting worth parameters
wmat <- llbt.worth(res)
plot(wmat)

# object specific covariates
LAT  <- c(0, 1, 1, 0, 1, 0)        # latin cities
EC   <- c(1, 0, 1, 0, 0, 1)
OBJ  <- data.frame(LAT, EC)
des2 <- llbt.design(cemspc, nitems = 6, objcovs = OBJ, cat.scovs = "ENG")
res2 <- gnm(y ~ LAT + EC + LAT:ENG + g1,
    eliminate = mu:ENG, data = des2, family = poisson)

# calculating and plotting worth parameters
wmat2 <- llbt.worth(res2)
wmat2
plot(wmat2)
}
\keyword{models}
\keyword{regression}
