## high-level convenience interface
raschtree <- function(formula, data, minsplit = 10, ...)
{
  ## transform formula
  stopifnot(length(formula) > 2)
  formula <-  formula(terms(formula, data = data))
  ff <- y ~ 1 | x
  ff[[2]] <- formula[[2]]
  ff[[3]][[3]] <- formula[[3]]

  ## formula/data/model pre-processing
  raschmod <- RaschModel()  
  ff <- attr(modeltools:::ParseFormula(ff), "formula")
  ff$input[[3]] <- ff$input[[2]]
  ff$input[[2]] <- ff$response[[2]]
  ff <- dpp(raschmod, as.formula(ff$input), other = list(part = as.formula(ff$blocks)), 
    data = data, na.action = na.pass)

  ## data sanity checks
  y <- as.matrix(ff@get("response"))
  if(ncol(y) < 3) stop("need at least three items")
  if(!all(as.vector(y) %in% c(0:1, NA))) stop("y must be a binary 0/1 matrix (potentially with NAs)")
  if(!all(apply(y, 1, function(x) all(0:1 %in% x))))
    stop("each row of y must have at least one 0 and one 1 entry")

  ## call mob()
  rval <- mob(ff, model = raschmod, control = mob_control(minsplit = minsplit,
      objfun = function(object) - as.vector(logLik(object)), ...))

  ## add class and return
  structure(list(mob = rval), class = "raschtree")
}

## convenience plotting
plot.raschtree <- function(x, terminal_panel = node_raschplot, tnex = 2, ...) {
  plot(x$mob, terminal_panel = terminal_panel, tnex = tnex, tp_args = list(...))
}

## hand-crafted "Next()" to bridge to
## un-exported S4 classes "mob"/"BinaryTree", argh!
logLik.raschtree <- function(object, ...) logLik(object$mob, ...)
sctest.raschtree <- function(x, ...) sctest(x$mob, ...)
weights.raschtree <- function(object, ...) weights(object$mob, ...)
summary.raschtree <- function(object, ...) summary(object$mob, ...)
print.raschtree <- function(x, ...) {
  print(x$mob, ...)
  invisible(x)
}

## parameters for Rasch trees
coef.raschtree <- function (object, node = NULL, ...) 
{
  object <- object$mob
  if(is.null(node)) node <- party:::terminal_nodeIDs(object@tree)
  rval <- sapply(nodes(object, node), function(z) coef(z$model, ...))
  if (!is.null(dim(rval))) {
    rval <- t(rval)
    rownames(rval) <- node
  }
  return(rval)
}

worth.raschtree <- function (object, node = NULL, ...) 
{
  object <- object$mob
  if(is.null(node)) node <- party:::terminal_nodeIDs(object@tree)
  rval <- sapply(nodes(object, node), function(z) worth(z$model, ...))
  if (!is.null(dim(rval))) {
    rval <- t(rval)
    rownames(rval) <- node
  }
  return(rval)
}

## visualization function
node_raschplot <- function(mobobj, id = TRUE,
  center = TRUE, names = TRUE, abbreviate = TRUE, index = TRUE, ref = TRUE,
  col = "black", linecol = "lightgray", cex = 0.5, pch = 19, xscale = NULL, yscale = NULL, ylines = 1.5)
{
    ## extract parameter of interest
    node <- 1:max(party:::terminal_nodeIDs(mobobj@tree))
    cf <- t(sapply(nodes(mobobj, node), function(z)
      if(center) worth(z$model) else worth(z$model) - worth(z$model)[1]))
    rownames(cf) <- node

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
    
        ## dependent variable setup
	cfi <- cf[node$nodeID,]
	cf_ref <- mean(cfi)

        ## viewport setup
        top_vp <- viewport(layout = grid.layout(nrow = 2, ncol = 3,
                           widths = unit(c(ylines, 1, 1), c("lines", "null", "lines")),  
			   heights = unit(c(1, 1), c("lines", "null"))),
                           width = unit(1, "npc"), 
                           height = unit(1, "npc") - unit(2, "lines"),
			   name = paste("node_raschplot", node$nodeID, sep = ""))
        pushViewport(top_vp)
        grid.rect(gp = gpar(fill = "white", col = 0))

        ## main title
        top <- viewport(layout.pos.col = 2, layout.pos.row = 1)
        pushViewport(top)
	mainlab <- paste(ifelse(id, paste("Node", node$nodeID, "(n = "), ""),
	                 sum(node$weights), ifelse(id, ")", ""), sep = "")
        grid.text(mainlab)
        popViewport()

        ## actual plot	
        plot_vpi <- viewport(layout.pos.col = 2, layout.pos.row = 2,
	    xscale = xscale, yscale = yscale, 
	    name = paste("node_raschplot", node$nodeID, "plot", sep = ""))
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

