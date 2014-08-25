## high-level convenience interface to mob()
raschtree <- function(formula, data, na.action = na.pass,
  reltol = 1e-10, deriv = c("sum", "diff", "numeric"), maxit = 100L, ...)
{
  ## keep call
  cl <- match.call(expand.dots = TRUE)

  ## use dots for setting up mob_control
  control <- mob_control(...)
  control$ytype <- "matrix"

  ## control options for raschfit
  raschcontrol <- list(reltol = reltol, deriv = deriv, maxit = maxit)

  ## call mob
  m <- match.call(expand.dots = FALSE)
  m$fit <- raschfit
  m$control <- control
  for(n in names(raschcontrol)) if(!is.null(raschcontrol[[n]])) m[[n]] <- raschcontrol[[n]]
  if("..." %in% names(m)) m[["..."]] <- NULL
  m[[1L]] <- as.name("mob")
  rval <- eval(m, parent.frame())
  
  ## extend class and keep original call
  rval$info$call <- cl
  class(rval) <- c("raschtree", class(rval))
  return(rval)
}

## glue code for calling RaschModel.fit()
raschfit <- function(y, x = NULL, start = NULL, weights = NULL, offset = NULL, ...,
  estfun = FALSE, object = FALSE)
{
  if(!(is.null(x) || NCOL(x) == 0L)) warning("x not used")
  if(!is.null(offset)) warning("offset not used")
  ## rval <- RaschModel.fit(y, weights = weights, start = start, ..., full = object | estfun)
  rval <- raschmodel(y, weights = weights, start = start, ..., hessian = object | estfun)
  ## rval <- RaschModel.fit(y, weights = weights, start = NULL, ..., full = object | estfun)
  ## rval <- RaschModel.fit(y, weights = weights, start = start, ..., full = TRUE)
  rval <- list(
    coefficients = rval$coefficients,
    objfun = -rval$loglik,
    estfun = if(estfun) estfun.raschmodel(rval) else NULL,
    object = if(object) rval else NULL
  )
  return(rval)
}

## methods
print.raschtree <- function(x,
  title = "Rasch tree", objfun = "negative log-likelihood", ...)
{
  partykit::print.modelparty(x, title = title, objfun = objfun, ...)
}

predict.raschtree <- function(object, newdata = NULL,
  type = c("worth", "rank", "best", "node"), ...)
{
  ## type of prediction
  type <- match.arg(type)
  
  ## nodes can be handled directly
  if(type == "node") return(partykit::predict.modelparty(object, newdata = newdata, type = "node", ...))
  
  ## get default newdata otherwise
  if(is.null(newdata)) newdata <- model.frame(object)
  
  pred <- switch(type,
    "worth" = worth,
    "rank" = function(obj, ...) rank(-worth(obj)),
    "best" = function(obj, ...) {
      wrth <- worth(obj)
      factor(names(wrth)[which.max(wrth)], levels = names(wrth))
    }
  )
  partykit::predict.modelparty(object, newdata = newdata, type = pred, ...)
}

apply_to_models <- function(object, node = NULL, FUN = NULL, drop = FALSE, ...) {
  if(is.null(node)) node <- nodeids(object, terminal = FALSE)
  if(is.null(FUN)) FUN <- function(object, ...) object  
  rval <- if("object" %in% object$info$control$terminal) {
    nodeapply(object, node, function(n) FUN(info_node(n)$object))
  } else {
    lapply(refit.modelparty(object, node, drop = FALSE), FUN)
  }
  names(rval) <- node
  if(drop & length(node) == 1L) rval <- rval[[1L]]
  return(rval)
}

itempar.raschtree <- function(object, node = NULL, ...)
{
  ids <- if(is.null(node)) nodeids(object, terminal = TRUE) else node
  if(length(ids) == 1L) {
    apply_to_models(object, node = ids, FUN = itempar, drop = TRUE)
  } else {
    do.call("rbind", apply_to_models(object, node = ids, FUN = itempar, drop = FALSE))
  } 
}

plot.raschtree <- function(x, type = "profiles", terminal_panel = NULL,
  tp_args = list(), tnex = NULL, drop_terminal = NULL, ...)
{
  if(!is.null(terminal_panel) && !missing(type)) {
    warning("Only one of 'type' and 'terminal_panel' should be specified")
  } else {
    terminal_panel <- switch(match.arg(type),
      "regions" = node_regionplot,
      "profiles" = node_profileplot)
  }
  if(is.null(tnex)) tnex <- if(is.null(terminal_panel)) 1L else 2L
  if(is.null(drop_terminal)) drop_terminal <- !is.null(terminal_panel)
  partykit::plot.modelparty(x, terminal_panel = terminal_panel,
    tp_args = tp_args, tnex = tnex, drop_terminal = drop_terminal, ...)
}
