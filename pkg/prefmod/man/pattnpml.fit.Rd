\name{pattnpml.fit}
\alias{pattnpml.fit}
\title{ NPML estimation for paired comparison models}
\description{ Fits a mixture model to overdispersed paired comparison data using nonparametric maximum likelihood (Aitkin, 1996a).
}
\usage{
pattnpml.fit(formula,
        random = ~1,
        k = 1,
        design,
        tol = 0.5,
        startp = NULL,
        EMmaxit = 500,
        EMdev.change = 0.001,
        pr.it = FALSE
        )
}
\arguments{
  \item{formula}{A formula defining the response (the count of the number of cases of each pattern) and the fixed effects (e.g. \code{y ~ x}).   }
  \item{random}{A formula defining the random model.  If there are three objects labelled o1, o2, o3, set \code{random = ~o1+o2+o3} to model overdispersion. For more details, see below.
  }
  \item{k}{The number of mass points (latent classes). Up to 21 mass points are supported.}
  \item{design}{The design data frame for paired comparison data as generated using  \code{\link{patt.design}} (mandatory, even if it is attached to the workspace!). }
  \item{tol}{The \code{tol} scalar (usually, \eqn{0 < }\code{tol} \eqn{\le 1}{<= 1}). This scalar sets the scaling factor for the locations of the initial mass points. A larger value means that the starting point locations are more widely spread}
  \item{startp}{Optional numerical vector of length \code{k} specifying the starting probabilities for the mass points to initiailise the EM algorithm. The default is to take gausssian quandrature probabilities  }
  \item{EMmaxit}{The maximum number of EM iterations.   }
  \item{EMdev.change}{Stops EM algorithm when deviance change falls below this value.  }
  \item{pr.it}{A dot is printed at each iteration cycle of the EM algorithm if set to \code{TRUE}.}
}
\details{
  The function \code{pattnpml.fit} is a wrapper function for \code{alldistPC} which
  in turn is a modified version of the function \code{alldist} from the
  \code{npmlreg} package.

  The nonparametric maximum likelihood (NPML) approach was introduced in Aitkin (1996)
  as a tool to fit overdispersed generalized linear models. The idea is to approximate
  the unknown and unspecified distribution of the random effect by a discrete mixture
  of exponential family densities, leading to a simple expression of the marginal
  likelihood which can then be maximized using a standard EM algorithm.

  This function extends the NPML aproach to allow fitting of overdispersed paired comparison models.
  It assumes that overdispersion arises because of dependence in the patterns. Fitting a non-parametric random effects
  term is equivalent to specifying distinct latent classes of response patterns.

  The number of components \code{k} of the finite mixture has to be specified beforehand.

  The EM algorithm used by the function takes the Gauss-Hermite masses
  and mass points as starting points. The position of the starting points can
  be concentrated or extended by setting \code{tol} smaller or larger,
  respectively; the initial mass point probabilities of the starting points can also be specified through \code{startp}.

  Fitting  models for overdispersion can be achieved by specifying the paired comparison items as additive terms
  in the random part of the model formula. A separate estimate for each item and for each mass point is produced.

  Fitting subject covariate models with the same effect for each mass
  point component is achieved by
  specifying as part of the \code{formula} a) a subject factor giving a
  different estimate for each covariate combination
  b) an interaction of the chosen subject covariates with the objects.  For
  models with subject factor covariates only, the first term
  is simply the interaction of all of the factor covariates.

  Fitting subject covariate models with a different effect for each mass
  point component (sometimes called random coefficient models, see Aitkin,
  Francis, Hinde and Darnell, 2009, pp. 497) is possible by specifying an
  interaction of the subject covariates with the items in the
  \code{random} term, and also in the \code{formula} part. Thus the setting \code{random= ~
  x:(o1+o2+o3} gives a model with a set of random slopes (one set for each
  mass point) and a set of random intercepts, one set for each mass point.

}
\value{
The function alldist produces an object of class \code{pattNPML}
The object contains the following 29 components:

  \item{coefficients}{a named vector of coefficients (including the mass points).
    In case of Gaussian quadrature, the coefficient given at \code{z}
    corresponds to the standard deviation of the mixing distribution. }
  \item{residuals}{the difference between the true response and the emprical Bayes predictions.}
  \item{fitted.values}{the empirical Bayes predictions (Aitkin, 1996b) on the scale of the responses.}
  \item{family}{the `family' object used.}
  \item{linear.predictors}{the extended linear predictors \eqn{\hat{\eta}_{ik}}{eta_ik}. }
  \item{disparity}{the disparity (\code{-2logL}) of the fitted mixture regression model.}
  \item{deviance}{the deviance of the fitted mixture regression model.}
  \item{null.deviance}{The deviance for the null model (just containing an intercept), comparable with
          `deviance'.}
  \item{df.residual}{the residual degrees of freedom of the fitted model (including the random part).}
  \item{df.null}{the residual degrees of freedom for the null model.}
  \item{y}{the (extended) response vector.}
  \item{call}{the matched call.}
  \item{formula}{the formula supplied.}
  \item{random}{the random term of the model formula.}
  \item{data}{the data argument.}
  \item{model}{the (extended) design matrix.}
  \item{weights}{the case weights initially supplied.}
  \item{offset}{the offset initially supplied.}
  \item{mass.points}{the fitted mass points.}
  \item{masses}{the mass point probabilities corresponding to the patterns.}
  \item{sdev}{a list of the two elements \code{sdev$sdev} and \code{sdev$sdevk}.
    The former is the estimated standard deviation of the
    Gaussian mixture components (estimated over all mixture components), and the latter
    gives the unequal or smooth component-specific standard deviations.
    All values are equal if \code{lambda=0}.}
  \item{shape}{a list of the two elements \code{shape$shape} and \code{shape$shapek},
    to be interpreted in analogy to sdev.}
  \item{rsdev}{estimated random effect standard deviation.}
  \item{post.prob}{a matrix of posteriori probabilities.}
  \item{post.int}{a vector of `posteriori intercepts' (as in Sofroniou et al. (2006)).}
  \item{ebp}{the empirical Bayes Predictions on the scale of the linear predictor.
  For compatibility with older versions.}
  \item{EMiter}{gives the number of iterations of the EM algorithm. }
  \item{EMconverged}{logical value indicating if the EM algorithm converged. }
  \item{lastglm}{the fitted \code{glm} object from the last EM iteration.}
  \item{Misc}{contains additional information relevant for the summary and plot functions,
   in particular the disparity trend and the EM trajectories. }

For further details see the help file for function \code{alldist} in package \code{npmlreg}.
}



\references{

Aitkin, M. (1996a). A general maximum likelihood analysis of overdispersion in generalized linear models. Statistics and Computing 6, 251-262.

Aitkin, M., Francis, B., Hinde, J. and Darnell, R. (2009). Statistical Modelling in R, Oxford Statistical Science Series, Oxford, UK.

Einbeck, J. & Hinde, J. (2006). A note on NPML estimation for exponential family regression models with unspecified dispersion parameter. Austrian Journal of Statistics 35, 233-243.

Sofroniou, N., Einbeck, J., and Hinde, J. (2006). Analyzing Irish suicide rates with mixture models. Proceedings of the 21st International Workshop on Statistical Modelling in Galway, Ireland, 2006.
}

\author{Originally translated from the GLIM 4 functions \code{alldist} and
\code{allvc} (Aitkin & Francis, 1995) to \R by Ross Darnell (2002). Modified,
extended, and prepared for publication by Jochen Einbeck and John Hinde (2006).
Adapted for paired comparison modelling by Reinhold Hatzinger and Brian Francis (2009)
}

\note{The mass point probabilities given in the output are the proportion of
patterns estimated to contribute to each mass point.
To estimate the  proportion  of cases contributing to each mass point the
posterior probabilities need to be averaged over patterns with observed
counts as weights (see example below).
}

\seealso{ \code{\link{glm}}.
}

\examples{
# two latent classes for paired comparison data
data(dat4)
dfr <- patt.design(dat4, 4)
modPC<-pattnpml.fit(y ~ 1, random = ~o1 + o2 + o3, k = 2, design = dfr)
modPC

# estimated proportion of cases in each mixture component
apply(modPC$post.prob, 2, function(x) sum(x * dfr$y / sum(dfr$y)))

\dontrun{
# fitting a model for two latent classes and fixed categorical
# subject covariates to the Eurobarometer 55.2 data (see help("euro55.2.des"))
# on rankings of sources of information on scientific developments

data(euro55.2.des)
model2cl <- pattnpml.fit(formula = y ~ -1 + SEX:AGE4 + (SEX + AGE4):
               (TV + RAD + NEWSP + SCIMAG + WWW + EDINST),
               random = ~ TV + RAD + NEWSP + SCIMAG + WWW + EDINST,
               k = 2, design = euro55.2.des, pr.it = T)
summary(model2cl)
}
}
\keyword{models}
\keyword{regression}
