if(FALSE) {
library("partykit")
library("psychotools")
data("Topmodel2007", package = "psychotree")
mb <-    mob(preference ~ age + gender + q1 + q2 + q3, data = Topmodel2007, fit = btfit)
bt <- bttree(preference ~ age + gender + q1 + q2 + q3, data = Topmodel2007)
}

## high-level convenience interface to mob()
bttree <- function(formula, data, na.action = na.pass,
  type = "loglin", ref = NULL, undecided = NULL, position = NULL, ...)
{
  ## keep call
  cl <- match.call(expand.dots = TRUE)

  ## use dots for setting up mob_control
  control <- mob_control(...)

  ## control options for btfit
  btcontrol <- list(type = type, ref = ref, undecided = undecided, position = position)

  ## call mob
  m <- match.call(expand.dots = FALSE)
  m$fit <- btfit
  m$control <- control
  for(n in names(btcontrol)) if(!is.null(btcontrol[[n]])) m[[n]] <- btcontrol[[n]]
  if("..." %in% names(m)) m[["..."]] <- NULL
  m[[1L]] <- as.name("mob")
  rval <- eval(m, parent.frame())
  
  ## extend class and keep original call
  rval$info$call <- cl
  class(rval) <- c("bttree", class(rval))
  return(rval)
}

## glue code for calling btReg.fit()
btfit <- function(y, x = NULL, start = NULL, weights = NULL, offset = NULL, ...,
  estfun = FALSE, object = FALSE)
{
  if(!(is.null(x) || NCOL(x) == 0L)) warning("x not used")
  if(!is.null(offset)) warning("offset not used")
  vcov <- object
  rval <- btReg.fit(y, weights = weights, start = start, ...,
    estfun = estfun, vcov = object)
  rval <- list(
    coefficients = rval$coefficients,
    objfun = -rval$loglik,
    estfun = rval$estfun,
    object = if(object) rval else NULL
  )
  return(rval)
}

## methods
print.bttree <- function(x,
  title = "Bradley-Terry tree", objfun = "negative log-likelihood", ...)
{
  partykit:::print.modelparty(x, title = title, objfun = objfun, ...)
}

predict.bttree <- function(object, newdata = NULL,
  type = c("worth", "rank", "best", "node"), ...)
{
  ## type of prediction
  type <- match.arg(type)
  
  ## nodes can be handled directly
  if(type == "node") return(partykit:::predict.modelparty(object, newdata = newdata, type = "node", ...))
  
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
  partykit:::predict.modelparty(object, newdata = newdata, type = pred, ...)
}

worth.bttree <- function(object, node = NULL, ...)
{
  ids <- if(is.null(node)) nodeids(object, terminal = TRUE) else node
  apply_to_models(object, node = ids, FUN = worth, drop = TRUE)
}

plot.bttree <- function(x, terminal_panel = node_btplot,
  tp_args = list(), tnex = NULL, drop_terminal = NULL, ...)
{
  if(is.null(tnex)) tnex <- if(is.null(terminal_panel)) 1L else 2L
  if(is.null(drop_terminal)) drop_terminal <- !is.null(terminal_panel)
  partykit:::plot.modelparty(x, terminal_panel = terminal_panel,
    tp_args = tp_args, tnex = tnex, drop_terminal = drop_terminal, ...)
}


## visualization function
node_btplot <- function(mobobj, id = TRUE,
  worth = TRUE, names = TRUE, abbreviate = TRUE, index = TRUE, ref = TRUE,
  col = "black", linecol = "lightgray", cex = 0.5, pch = 19, xscale = NULL, yscale = NULL, ylines = 1.5)
{
    ## node ids
    node <- nodeids(mobobj@tree, terminal = FALSE)
    
    ## get all coefficients 
    cf <- apply_to_models(mobobj, node, FUN = function(z)        
      if(worth) worth(z$model) else coef(z$model, all = FALSE, ref = TRUE))
    cf <- do.call("rbind", cf)
    rownames(cf) <- node

    ## get one full model
    mod <- apply_to_model(mobobj, node = 1L, FUN = NULL)

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
         
    ## panel function for bt plots in nodes
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
        		 name = paste("node_btplot", id, sep = ""))
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
        name = paste("node_btplot", id, "plot", sep = ""))
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
class(node_btplot) <- "grapcon_generator"
