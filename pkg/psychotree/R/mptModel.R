## S4 StatModel object
mptModel <- function(mptform = NULL, mptstruc = NULL, maxit = 1000) {
  new("StatModel",
    capabilities = new("StatModelCapabilities"),
    name = "Multinomial processing tree model",
    dpp = ModelEnvFormula,
    fit = function(object, weights = NULL, ...){

      ## extract response (there are no regressors)
      y <- object@get("response")[[1]]

      ## call btReg.fit()
      z <- mptModel.fit(y = y, weights = weights, mptform = mptform,
                        mptstruc = mptstruc, maxit = maxit)
      z$ModelEnv <- object
      z$addargs <- list(...)
      z
    }
  )
}


reweight.mptModel <- function(object, weights, ...) {
     fit <- mptModel(mptform = object$mptform, mptstruc = object$mptstruc,
                     maxit = object$maxit)@fit
     do.call("fit", c(list(object = object$ModelEnv, weights = weights),
                      object$addargs))
}


estfun.mptModel <- function(x, ...) x$estfun


# estfun.mptModel <- function(object, ...){
#   y       <- object$y
#   theta   <- coef(object)
#   pcat    <- object$pcat
#   sympcat <- parse(text = as.character(object$mptform[-1]))
# 
#   ## Symbolic derivative of sympcat for each theta
#   dpcat <- matrix(NA, nrow=length(sympcat), ncol=length(theta))
#   for(i in seq_along(sympcat))
#     for(j in seq_along(theta))
#       dpcat[i, j] <- eval(D(sympcat[i], names(theta)[j]), as.list(theta))
# 
#   # (Symbolic pcat evaluated at theta) == pcat
#   # pcat <- sapply(sympcat, eval, as.list(theta))
# 
#   # Wide data
#   t(sapply(seq_len(nrow(y)), function(i) colSums(y[i,]/pcat * dpcat)))
# }


bread.mptModel <- function(x, ...) vcov(x) * x$n

