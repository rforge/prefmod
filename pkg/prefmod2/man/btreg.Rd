\name{btreg}

\alias{btreg}
\alias{btreg.fit}

\alias{print.btreg}
\alias{summary.btreg}
\alias{print.summary.btreg}

\alias{coef.btreg}
\alias{worth.btreg}
\alias{deviance.btreg}
\alias{logLik.btreg}
\alias{vcov.btreg}
\alias{bread.btreg}
\alias{estfun.btreg}

\title{Bradley-Terry Regression Models}

\description{
  Extended Bradley-Terry models with subject- and object-specific
  variables fitted by means of log-linear (or logistic) regression
  models.
}

\usage{
btreg(formula, data, subset, na.action, weights, offset,
  type = c("loglin", "logit"), ref = NULL, undecided = NULL, position = NULL,
  model = FALSE, y = TRUE, x = FALSE, \dots)

btreg.fit(x, y, z = NULL, weights = NULL, type = c("loglin", "logit"), ref = NULL,
  undecided = NULL, position = NULL, \dots)
}

\arguments{
  \item{formula}{A symbolic description of the model to be fit. This
    should be of type \code{y | z1 + z2 ~ x1 + x2 + x3}
    where \code{y} should be an object of class \code{\link[psychotree]{paircomp}},
    \code{x1}, \code{x2}, \code{x3} are subject-covariates (taken from \code{data})
    and \code{z1} and \code{z2} are object-covariates (taken from \code{covariates(y)}).}
  \item{data}{an optional data frame containing the variables in the model.}
  \item{subset}{an optional subset specification.}
  \item{na.action}{A function which indicates what should happen when the data
    contain \code{NA}s, defaulting to \code{\link[stats]{na.action}}.}
  \item{weights}{an optional vector of weights, interpreted as case weights (integer only).}
  \item{offset}{an optional offset vector (currently ignored).}
  \item{type}{character. Should an auxiliary log-linear Poisson model or logistic
    binomial be employed for estimation? The latter is only available if not undecided
    effects are estimated.}
  \item{ref}{character or numeric. Which object parameter should be the reference
    category, i.e., constrained to zero?}
  \item{undecided}{logical. Should an undecided parameter be estimated?}
  \item{position}{logical. Should a position effect be estimated?}
  \item{model}{logical. Should the model frame be included in the
    fitted model object?}
  \item{x}{matrix. The regressor matrix pertaining to the subject covariates.
    This can be returned in \code{btreg} if \code{x = TRUE}.}
  \item{y}{paircomp object with the response.
    This can be returned in \code{btreg} if \code{y = TRUE}.}
  \item{z}{matrix. The regressor matrix pertaining to the object covariates.
    This can be returned in \code{btreg} if \code{x = TRUE}.}
  \item{\dots}{further arguments passed to functions.}
}

\details{
  \code{btreg} is a convenient formula-based interface for fitting Bradley-Terry
  models. It essentially handles the formula parsing and argument preprocessing
  and then calls \code{btreg.fit}, the workhorse fitting function. This suitably
  aggregates the data and then calls \code{\link{glm.fit}} (currently) for parameter
  estimation and postprocesses the results.
  
  Currently, the fitting function only supports simple Bradley-Terry models
  of type \code{y ~ 1}, i.e., without any covariates. This will be extended in
  future versions.
}

\value{
  An object of class \code{"btreg"}, i.e., a list with components including
  \item{coefficients}{a vector of estimated coefficients,}
  \item{vcov}{covariance matrix of all coefficients in the model,}
  \item{loglik}{log-likelihood of the fitted model,}
  \item{df}{degrees of freedom,}  
  \item{estfun}{empirical estimating functions (gradients),}
  \item{weights}{case weights used,}
  \item{n}{number of observations (subjects),}
  \item{type}{\code{type} argument used,}
  \item{ref}{\code{ref} argument used,}
  \item{undecided}{\code{undecided} argument used,}
  \item{position}{\code{position} argument used,}
  \item{labels}{vector of object labels,}
  \item{call}{the original function call,}
  \item{terms}{the original model terms,}
  \item{model}{the full model frame (if \code{model = TRUE}),}
  \item{y}{the response paircomp vector (if \code{y = TRUE}),}
  \item{x}{regressor matrix of subject covariates (if \code{x = TRUE}),}
  \item{z}{regressor matrix of object covariates (if \code{x = TRUE}),}
}

\references{ 
  Critchlow DE, Fligner MA (1991). Paired Comparison, Triple Comparison, and Ranking
    Experiments as Generalized Linear Models, and their Implementation in GLIM.
    \emph{Psychometrika}, \bold{56}(3), 517--533.

  Dittrich R, Hatzinger R, Katzenbeisser W (1998). Modelling the Effect of
    Subject-Specific Covariates in Paired Comparison Studies with an Application
    to University Rankings.
    \emph{Journal of the Royal Statistical Society C}, \bold{47}(4), 511--525.
}

\seealso{\code{\link{plot.btreg}}}

\examples{
## CEMS university choice data
data("CEMSChoice", package = "psychotree")

## Bradley-Terry modeling
cems <- btreg(preference ~ 1, data = CEMSChoice)
cems
summary(cems)
plot(cems)
worth(cems)
}

\keyword{regression}
