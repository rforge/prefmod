## high-level convenience interface to mob()
rstree <- function(formula, data, na.action = na.pass,
  reltol = 1e-10, deriv = c("sum", "diff"), maxit = 100L, ...)
{
  ## keep call
  cl <- match.call(expand.dots = TRUE)

  ## use dots for setting up mob_control
  control <- mob_control(...)
  control$ytype <- "matrix"

  ## control options for rsmfit
  rsmcontrol <- list(reltol = reltol, deriv = deriv, maxit = maxit)

  ## call mob
  m <- match.call(expand.dots = FALSE)
  m$fit <- rsmfit
  m$control <- control
  for(n in names(rsmcontrol)) if(!is.null(rsmcontrol[[n]])) m[[n]] <- rsmcontrol[[n]]
  if("..." %in% names(m)) m[["..."]] <- NULL
  m[[1L]] <- as.name("mob")
  rval <- eval(m, parent.frame())
  
  ## extend class and keep original call
  rval$info$call <- cl
  class(rval) <- c("rstree", class(rval))
  return(rval)
}

## glue code for calling RSModel.fit()
rsmfit <- function(y, x = NULL, start = NULL, weights = NULL, offset = NULL, ...,
  estfun = FALSE, object = FALSE)
{
  if(!(is.null(x) || NCOL(x) == 0L)) warning("x not used")
  if(!is.null(offset)) warning("offset not used")
  rval <- rsmodel(y, weights = weights, start = start, ..., hessian = object | estfun)
  rval <- list(
    coefficients = rval$coefficients,
    objfun = -rval$loglik,
    estfun = if(estfun) estfun.rsmodel(rval) else NULL,
    object = if(object) rval else NULL
  )
  return(rval)
}

## methods
print.rstree <- function(x,
  title = "Rating scale tree", objfun = "negative log-likelihood", ...)
{
  partykit::print.modelparty(x, title = title, objfun = objfun, ...)
}

predict.rstree <- function(object, newdata = NULL,
  type = c("node"), ...)
{
  ## type of prediction
  type <- match.arg(type)
  
  ## nodes can be handled directly
  if(type == "node") return(partykit::predict.modelparty(object, newdata = newdata, type = "node", ...))
  
  ## ## get default newdata otherwise
  ## if(is.null(newdata)) newdata <- model.frame(object)
  
  ## pred <- switch(type,
  ##   "worth" = worth,
  ##   "rank" = function(obj, ...) rank(-worth(obj)),
  ##   "best" = function(obj, ...) {
  ##     wrth <- worth(obj)
  ##     factor(names(wrth)[which.max(wrth)], levels = names(wrth))
  ##   }
  ## )
  ## partykit::predict.modelparty(object, newdata = newdata, type = pred, ...)
}

## worth.rstree <- function(object, node = NULL, ...)
## {
##   ids <- if(is.null(node)) nodeids(object, terminal = TRUE) else node
##   if(length(ids) == 1L) {
##     apply_to_models(object, node = ids, FUN = worth, drop = TRUE)
##   } else {
##     do.call("rbind", apply_to_models(object, node = ids, FUN = worth, drop = FALSE))
##   } 
## }

plot.rstree <- function(x, terminal_panel = node_regionplot,
  tp_args = list(), tnex = NULL, drop_terminal = NULL, ...)
{
  if(is.null(tnex)) tnex <- if(is.null(terminal_panel)) 1L else 2L
  if(is.null(drop_terminal)) drop_terminal <- !is.null(terminal_panel)
  partykit::plot.modelparty(x, terminal_panel = terminal_panel,
    tp_args = tp_args, tnex = tnex, drop_terminal = drop_terminal, ...)
}
