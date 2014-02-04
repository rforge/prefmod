if(FALSE) {
library("partykit")
library("psychotools")
data("DIFSim", package = "psychotree")
mb <-       mob(resp ~ age + gender + motivation, data = DIFSim, fit = raschfit, control = mob_control(ytype = "matrix"))
rt <- raschtree(resp ~ age + gender + motivation, data = DIFSim)
}

## high-level convenience interface to mob()
raschtree <- function(formula, data, na.action = na.pass,
  gradtol = 1e-6, deriv = c("sum", "diff", "numeric"), ...)
{
  ## keep call
  cl <- match.call(expand.dots = TRUE)

  ## use dots for setting up mob_control
  control <- mob_control(...)
  control$ytype <- "matrix"

  ## control options for raschfit
  raschcontrol <- list(gradtol = gradtol, deriv = deriv)

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
  rval <- RaschModel.fit(y, weights = weights, start = start, ..., hessian = object | estfun)
  ## rval <- RaschModel.fit(y, weights = weights, start = NULL, ..., full = object | estfun)
  ## rval <- RaschModel.fit(y, weights = weights, start = start, ..., full = TRUE)
  rval <- list(
    coefficients = rval$coefficients,
    objfun = -rval$loglik,
    estfun = if(estfun) estfun.RaschModel(rval) else NULL,
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

worth.raschtree <- function(object, node = NULL, ...)
{
  ids <- if(is.null(node)) nodeids(object, terminal = TRUE) else node
  if(length(ids) == 1L) {
    apply_to_models(object, node = ids, FUN = worth, drop = TRUE)
  } else {
    do.call("rbind", apply_to_models(object, node = ids, FUN = worth, drop = FALSE))
  } 
}

plot.raschtree <- function(x, terminal_panel = node_raschplot,
  tp_args = list(), tnex = NULL, drop_terminal = NULL, ...)
{
  if(is.null(tnex)) tnex <- if(is.null(terminal_panel)) 1L else 2L
  if(is.null(drop_terminal)) drop_terminal <- !is.null(terminal_panel)
  partykit::plot.modelparty(x, terminal_panel = terminal_panel,
    tp_args = tp_args, tnex = tnex, drop_terminal = drop_terminal, ...)
}


## visualization function
node_raschplot <- function(mobobj, id = TRUE,
  worth = TRUE, names = TRUE, abbreviate = TRUE, index = TRUE, ref = TRUE,
  col = "black", linecol = "lightgray", cex = 0.5, pch = 19, xscale = NULL, yscale = NULL, ylines = 1.5)
{
    ## node ids
    node <- nodeids(mobobj, terminal = FALSE)
    
    ## get all coefficients 
    cf <- apply_to_models(mobobj, node, FUN = function(z)        
      if(worth) worth(z) else coef(z, all = FALSE, ref = TRUE))
    cf <- do.call("rbind", cf)
    rownames(cf) <- node

    ## get one full model
    mod <- apply_to_models(mobobj, node = 1L, FUN = NULL)

    if(!worth) {
      if(is.character(ref) | is.numeric(ref)) {
        reflab <- ref
        ref <- TRUE
      } else {
        reflab <- mod$ref
      }
      if(is.character(reflab)) reflab <- match(reflab, mod$labels)
      cf <- cf - cf[,reflab]
    }

    ## reference
    if(worth) {
      cf_ref <- 1/ncol(cf)
    } else {
      cf_ref <- 0
    }

    ## labeling
    if(is.character(names)) {
      colnames(cf) <- names
      names <- TRUE
    }

    ## abbreviation
    if(is.logical(abbreviate)) {
      nlab <- max(nchar(colnames(cf)))
      abbreviate <- if(abbreviate) as.numeric(cut(nlab, c(-Inf, 1.5, 4.5, 7.5, Inf))) else nlab
    }
    colnames(cf) <- abbreviate(colnames(cf), abbreviate)
    
    if(index) {
      x <- 1:NCOL(cf)
      if(is.null(xscale)) xscale <- range(x) + c(-0.1, 0.1) * diff(range(x))
    } else {
      x <- rep(0, length(cf))
      if(is.null(xscale)) xscale <- c(-1, 1)      
    }
    if(is.null(yscale)) yscale <- range(cf) + c(-0.1, 0.1) * diff(range(cf))
         
    ## panel function for rasch plots in nodes
    rval <- function(node) {

      ## node index
      id <- id_node(node)
    
      ## dependent variable setup
      cfi <- cf[id,]

      ## viewport setup
      top_vp <- viewport(layout = grid.layout(nrow = 2, ncol = 3,
    			 widths = unit(c(ylines, 1, 1), c("lines", "null", "lines")),  
        		 heights = unit(c(1, 1), c("lines", "null"))),
    			 width = unit(1, "npc"), 
    			 height = unit(1, "npc") - unit(2, "lines"),
        		 name = paste("node_raschplot", id, sep = ""))
      pushViewport(top_vp)
      grid.rect(gp = gpar(fill = "white", col = 0))

      ## main title
      top <- viewport(layout.pos.col = 2, layout.pos.row = 1)
      pushViewport(top)
      mainlab <- paste(ifelse(id, paste("Node", id, "(n = "), ""),
        	       info_node(node)$nobs, ifelse(id, ")", ""), sep = "")
      grid.text(mainlab)
      popViewport()

      ## actual plot  
      plot_vpi <- viewport(layout.pos.col = 2, layout.pos.row = 2,
        xscale = xscale, yscale = yscale, 
        name = paste("node_raschplot", id, "plot", sep = ""))
      pushViewport(plot_vpi)

      grid.lines(xscale, c(cf_ref, cf_ref), gp = gpar(col = linecol), default.units = "native")
      if(index) {
        grid.lines(x, cfi, gp = gpar(col = col, lty = 2), default.units = "native")
        grid.points(x, cfi, gp = gpar(col = col, cex = cex), pch = pch, default.units = "native")
        grid.xaxis(at = x, label = if(names) names(cfi) else x)
      } else {  	
        if(names) grid.text(names(cfi), x = x, y = cfi, default.units = "native")
          else grid.points(x, cfi, gp = gpar(col = col, cex = cex), pch = pch, default.units = "native")
      }
      grid.yaxis(at = c(ceiling(yscale[1] * 100)/100, floor(yscale[2] * 100)/100))
      grid.rect(gp = gpar(fill = "transparent"))

      upViewport(2)
    }
	    
    return(rval)
}
class(node_raschplot) <- "grapcon_generator"


bread.RaschModel <- function(x, ...) x$vcov * x$n

estfun.RaschModel <- function(x, ...) {
  ## extract data and parameters of interest
  par <- x$coefficients
  esf <- x$elementary_symmetric_functions
  y <- x$data
  weights_orig <- weights(x)
  y <- y[weights_orig > 0, , drop = FALSE]
  weights <- weights_orig[weights_orig > 0]
  rs <- rowSums(y)
  
  ## analytical gradient
  if(!x$na) {
    agrad <- weights * (- y + esf[[2]][rs + 1, , drop = FALSE] / esf[[1]][rs + 1])[,-1, drop = FALSE]
  } else {
    ## set up return value
    n <- nrow(y)
    k <- ncol(y)
    agrad <- matrix(0, nrow = n, ncol = k)

    ## observed NA patterns
    na_patterns <- factor(apply(is.na(y), 1, function(z) paste(which(z), collapse = "\r")))

    ## loop over observed NA patterns	   
    for(i in seq_along(levels(na_patterns))) {
      ## parse NA pattern
      lev_i <- levels(na_patterns)[i]
      na_i <- which(na_patterns == lev_i)
      wi_i <- as.integer(strsplit(lev_i, "\r")[[1]])
      wi_i <- if(length(wi_i) < 1) 1:k else (1:k)[-wi_i]

      ## compute gradient per pattern
      esf_i <- esf[[i]]
      rs_i <- rowSums(y[na_i, wi_i, drop = FALSE])
      agrad[na_i, wi_i] <- weights[na_i] * (- y[na_i, wi_i, drop = FALSE] +
    	esf_i[[2]][rs_i + 1, , drop = FALSE] / esf_i[[1]][rs_i + 1])
    }

    agrad <- agrad[, -1, drop = FALSE]
  }

  ## collect and return
  grad <- matrix(0, ncol = length(par), nrow = length(weights_orig))
  grad[weights_orig > 0,] <- agrad
  return(grad)
}

