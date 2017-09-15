\encoding{UTF-8}
\name{llbt.fit}
\alias{llbt.fit}
\title{Function to fit an LLBT}
\description{%
  Function to fit an LLBT using an ELIMINATE feature
}
\usage{
llbt.fit(y, Xmodel, q, ncat, maxiter = 100)
}
\arguments{%
  \item{y}{response, usually counts}
  \item{Xmodel}{design matrix}
  \item{q}{number of parameters to eliminate (usually number of comparisons times number of subject covariate levels}
  \item{ncat}{number of response categories}
  \item{maxiter}{maximum number of iterations (default 100)}
}
\details{%
  Be careful when specifying the design matrix.
  Since there is no extrinsic aliasing the matrix must have full rank.
  Usually, one of the design columns for object must be left out.
}
\references{%
  Hatzinger, R., & Francis, B. (2004). \emph{Fitting paired comparison models in R.} \url{http://epub.wu.ac.at/id/eprint/740}
}
\author{Reinhold Hatzinger}
\examples{
# fit basic model casewise
mfr <- llbt.design(cemspc, nitems = 6,
    objnames = c("lo", "pa", "mi", "sg", "ba", "st"),
    casewise = TRUE)
mm   <- model.matrix(~ lo+pa+mi+sg+ba + g1, data = mfr)
X    <- mm[, -1]
p    <- ncol(X)
ncat <- 3
q    <- length(levels(mfr$mu)) * length(levels(mfr$CASE))
llbt.fit(mfr$y, X, q, ncat)


# fit the (aggregated) model with one subject covariate
mfr <- llbt.design(cemspc, nitems = 6,
    objnames = c("lo", "pa", "mi", "sg", "ba", "st"),
    cov.sel = "ENG")
eng  <- mfr$ENG
eng  <- factor(eng)
mm   <- model.matrix(~ lo+pa+mi+sg+ba + g1 + (lo+pa+mi+sg+ba):eng, data = mfr)
X    <- mm[, -1]
q    <- length(levels(mfr$mu)) * length(levels(eng))
ncat <- 3
llbt.fit(mfr$y, X, q, ncat)
}
\keyword{models}
