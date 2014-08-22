##### panel generating graphic functions for various tree models

## wrapper functions for backward compatibility
node_raschplot <- function (...)
{
  node_profileplot(...)
}
class(node_raschplot) <- "grapcon_generator"


node_effects <- function (...)
{
  node_regionplot(...)
}
class(node_effects) <- "grapcon_generator"


## profile plot visualization function
node_profileplot <- function(mobobj, what = c("items", "thresholds", "discriminations"),
  paramarg = list(type = NULL, ref = NULL, alias = TRUE), id = TRUE, names = TRUE,
  abbreviate = TRUE, index = TRUE, ref = TRUE, col = "black", linecol = "lightgray",
  cex = 0.5, pch = 19, xscale = NULL, yscale = NULL, ylines = 1.5, ...)
{
  ## check input
  what <- match.arg(what)
  if (what == "thresholds") type <- paramarg$type
  refpar <- paramarg$ref
  alias <- if (is.null(paramarg$alias)) TRUE else paramarg$alias
  addargs <- list(...)
  if ("worth" %in% names(addargs)) warning("The argument 'worth' is deprecated and not longer used.")

  ## node ids
  node <- nodeids(mobobj, terminal = FALSE)
  
  ## get all coefficients 
  if (what == "items") {
    cf <- apply_to_models(mobobj, node, FUN = function(z) coef(itempar(z, ref = refpar, alias = alias, vcov = FALSE)))
  } else if (what == "thresholds") {
    cf <- apply_to_models(mobobj, node, FUN = function(z) coef(threshpar(z, type = type, ref = refpar, alias = alias, vcov = FALSE), type = "matrix"))
  } else {
    cf <- apply_to_models(mobobj, node, FUN = function(z) coef(discrpar(z, ref = refpar, alias = alias, vcov = FALSE)))
  }
  names(cf) <- node

  ## labeling
  if (is.logical(names) & names) {
    if (what != "thresholds") nms <- lapply(cf, names) else nms <- lapply(cf, function (j) rownames(j))
  } else if (is.character(names)) {
    nms <- split(rep(names, length(node)), f = rep(1:length(node), each = length(names)))
    names <- TRUE
  } else {
    names <- FALSE
  }

  ## abbreviation
  if (is.logical(abbreviate)) {
    nlab <- max(unlist(lapply(nms, function (j) nchar(j))))
    abbreviate <- if (abbreviate) as.numeric(cut(nlab, c(-Inf, 1.5, 4.5, 7.5, Inf))) else nlab
  }
  nms <- lapply(nms, function (j) abbreviate(j, abbreviate))

  ## axis scale
  if (index) {
    x <- if (what == "thresholds") 1:nrow(cf[[1]]) else 1:length(cf[[1]])
    if (is.null(xscale)) xscale <- range(x) + c(-0.1, 0.1) * diff(range(x))
  } else {
    if (what == "thresholds") {
      x <- 0:(ncol(cf[[1]]) - 1)
      if (is.null(xscale)) xscale <- c(-1, ncol(cf[[1]]))
    } else {
      x <- rep(0, length(cf[[1]]))
      if (is.null(xscale)) xscale <- c(-1, 1)
    }
  }
  r <- diff(range(unlist(cf), na.rm = TRUE))
  if (!r) r <- 1
  if (is.null(yscale)) yscale <- range(unlist(cf), na.rm = TRUE) + c(-0.1, 0.1) * r

  ## panel function for profile plots in nodes
  panelfun <- function (node) {

    ## node index
    idn <- id_node(node)
    
    ## get cfs and labels
    cfi <- cf[[idn]]
    if (names) nmsi <- nms[[idn]]


    ## viewport setup
    top_vp <- grid::viewport(layout = grid::grid.layout(nrow = 2, ncol = 3,
                           widths = grid::unit(c(ylines, 1, 1), c("lines", "null", "lines")),  
                           heights = grid::unit(c(1, 1), c("lines", "null"))),
                       width = grid::unit(1, "npc"), height = grid::unit(1, "npc") - grid::unit(2, "lines"),
                       name = paste("node_profileplot", idn, sep = ""))
    grid::pushViewport(top_vp)
    grid::grid.rect(gp = grid::gpar(fill = "white", col = 0))

    ## main title
    top <- grid::viewport(layout.pos.col = 2, layout.pos.row = 1)
    grid::pushViewport(top)
    mainlab <- paste(ifelse(id, paste("Node", idn, "(n = "), ""),
                     info_node(node)$nobs, ifelse(idn, ")", ""), sep = "")
    grid::grid.text(mainlab)
    grid::popViewport()

    ## actual plot  
    plot_vpi <- grid::viewport(layout.pos.col = 2, layout.pos.row = 2, xscale = xscale, yscale = yscale, 
                               name = paste("node_profileplot", idn, "plot", sep = ""))
    grid::pushViewport(plot_vpi)
    grid::grid.lines(xscale, c(mean(cfi), mean(cfi)), gp = grid::gpar(col = linecol), default.units = "native")
    if(index) {
      if (what == "thresholds") {
        for (j in 1:ncol(cfi)) {
          grid::grid.lines(x, cfi[, j], gp = grid::gpar(col = col, lty = 2), default.units = "native")
          grid::grid.text(label = paste0("C", j), x, cfi[, j], gp = grid::gpar(col = col), default.units = "native")
        }
      } else {
        grid::grid.lines(x, cfi, gp = grid::gpar(col = col, lty = 2), default.units = "native")
        grid::grid.points(x, cfi, gp = grid::gpar(col = col, cex = cex), pch = pch, default.units = "native")
      }
      grid::grid.xaxis(at = x, label = if (names) nmsi else x)
    } else {	
      if (names) {
        if (what == "thresholds") {
          for (j in 1:ncol(cfi)) grid::grid.text(paste0(nmsi, "-C", j), x[j], y = cfi[, j], default.units = "native")
        } else {
          grid::grid.text(nmsi, x = x, y = cfi, default.units = "native")
        }
      } else {
        if (what == "thresholds") {
          for (j in 1:ncol(cfi)) grid::grid.points(x = rep(j, nrow(cfi)), y = cfi[, j],
                                                   gp = grid::gpar(col = col, cex = cex), pch = pch, default.units = "native")
        } else {
          grid::grid.points(x, cfi, gp = grid::gpar(col = col, cex = cex), pch = pch, default.units = "native")
        }
      }
    }
    grid::grid.yaxis(at = c(ceiling(yscale[1] * 100)/100, floor(yscale[2] * 100)/100))
    grid::grid.rect(gp = grid::gpar(fill = "transparent"))

    grid::upViewport(2)
  }
  
  return(panelfun)
}
class(node_profileplot) <- "grapcon_generator"


## region plot visualization function
node_regionplot <- function(mobobj, names = NULL, type = c("mode", "median", "mean"),
  ref = NULL, ylab = "Latent trait", ylim = NULL, off = 0.1, col_fun = gray.colors,
  uo_show = TRUE, uo_col = "red", uo_lty = 2, uo_lwd = 1.25)    
{
    ## check input
    stopifnot(!is.null(mobobj))
    stopifnot(off >= 0)
    type <- match.arg(type)

    ## function to extract absolute item threshold parameters from model objects in terminal nodes
    threshparlst <- function (node) threshpar(node, ref = ref, type = type, alias = TRUE, relative = FALSE, cumulative = FALSE, vcov = FALSE)

    ## setup threshold parameters
    node <- nodeids(mobobj, terminal = TRUE)
    delta_lst <- apply_to_models(mobobj, node, FUN = threshparlst, ref = ref, type = type)

    ## setup plotting parameters
    m <- max(sapply(delta_lst, length))
    xi <- 0:m + c(0:(m - 1), m - 1) * off
    xlim <- c(xi[1], xi[m + 1])

    ## setup axis range and labels
    if (is.null(ylim)) ylim <- extendrange(unlist(delta_lst, use.names = FALSE), f = 0.25)
    if (is.null(names)) names  <- paste("I", 1:m, sep = "") else stopifnot(length(names) == m)

    ## label for extraction
    names(delta_lst) <- node

    ## finally, panel function, just select and paint.
    panelfun <- function(node) {

        ## select node id, coefficients, identified items, x-position vector and label
        id <- as.character(id_node(node))
        delta_unsorted <- delta_lst[[id]]
        lab <- paste("node", id, sep = "")
        
        ## compute sorted absolute item threshold parameters (code borrowed from psychotools::regionplot())
        delta_sorted <- delta_unsorted
        us <- sapply(delta_unsorted, is.unsorted)
        if (any(us)) {
          usj <- which(us)
          for (j in usj) {
            tpj <- delta_unsorted[[j]]
            nj <- length(tpj)
            
            ## check if there is a point with a biggest parameter, if yes, take mean
            for (i in 1:nj) {
              if (all(tpj[i] > tpj[(i+1):nj])) {
                tpj[i] <- mean(tpj[i:nj])
                tpj <- tpj[-(i+1:nj)]
                break
              }
            }
            
            ## recursive sorting if there is still unorder (e.g. 4, 2, 3, 1)
            while(is.unsorted(tpj)) {
              uo_pos <- which(diff(tpj) < 0)                     # locate unordered parameters, returns position of the first
              tpj[uo_pos] <- (tpj[uo_pos] + tpj[uo_pos + 1]) / 2 # replace first with location of intersection of ccc curves (= (eps1 + eps2)/ 2)
              tpj <- tpj[-(uo_pos + 1)]                          # remove second
            }
            
            delta_sorted[[j]] <- tpj
          }
        }

        ## terminal panel viewport setup
        top.vp <- grid::viewport(layout = grid::grid.layout(nrow = 2, ncol = 1, widths = grid::unit(1, "null"), heights = grid::unit(c(1, 1), c("lines", "null"))),
                                 width = grid::unit(1, "npc"), height = grid::unit(1, "npc") - grid::unit(2, "lines"), name = paste(lab, "_effects", sep = ""))
        grid::pushViewport(top.vp)
        grid::grid.rect(gp = grid::gpar(fill = "white", col = 0), name = paste(lab, "_border", sep = ""))

        ## main title
        grid::pushViewport(grid::viewport(layout.pos.col = 1, layout.pos.row = 1, name = paste(lab, "_title_vp", sep = "")))
        grid::grid.text(paste("Node ", id, " (n = ", info_node(node)$nobs, ")", sep = ""), name = paste(lab, "_title", sep = ""))
        grid::upViewport()
        
        ## finally the actual plot
        grid::pushViewport(grid::viewport(layout.pos.col = 1, layout.pos.row = 2, name = lab))
        lab <- paste(lab, "_plot", sep = "")

        ## setup plotting area (3x3)
        wcol <- if (is.null(ylab)) c(2.5, 1, 1) else c(4, 1, 1)
        hrow <- c(0.5, 1, 1)

        top.vp <- grid::viewport(layout = grid::grid.layout(nrow = 3, ncol = 3, widths = grid::unit(wcol, c("lines", "null", "lines")), heights = grid::unit(hrow, c("lines", "null", "lines"))), name = paste(lab, "_top_vp", sep = ""))
        bmargin.vp <- grid::viewport(layout.pos.row = 3, layout.pos.col = 2, name = paste(lab, "_bottom-margin_vp", sep = ""))
        lmargin.vp <- grid::viewport(layout.pos.row = 2, layout.pos.col = 1, name = paste(lab, "_left-margin_vp", sep = ""))
        rmargin.vp <- grid::viewport(layout.pos.row = 2, layout.pos.col = 3, name = paste(lab, "_right-margin_vp", sep = ""))
        plot.vp <- grid::viewport(layout.pos.row = 2, layout.pos.col = 2, name = paste(lab, "_vp", sep = ""), xscale = xlim, yscale = ylim)
        grid::pushViewport(top.vp)
        grid::pushViewport(plot.vp)

        ## plot rectangles per item
        for (j in seq_along(delta_sorted)) {
            ncat <- length(delta_sorted[[j]]) + 1
            grid::grid.rect(x = rep.int(xi[j], ncat), y = c(ylim[1], delta_sorted[[j]]), width = rep.int(1, ncat),
                            height = diff.default(c(ylim[1], delta_sorted[[j]], ylim[2])), just = c("left", "bottom"),
                            gp = grid::gpar(fill = col_fun(ncat)), default.units = "native", name = paste(lab, "_item", j, "_rect", sep = ""))
        }

        ## if requested: indicate unordered parameters
        if (uo_show && type == "mode") {
          uo_items <- which(!sapply(mapply(all.equal, delta_sorted, delta_unsorted, check.attributes = FALSE, SIMPLIFY = FALSE, USE.NAMES = FALSE), is.logical))
          for (j in uo_items) {
            uo_pars <- setdiff(delta_unsorted[[j]], delta_sorted[[j]])
            grid::grid.polyline(x = rep(c(xi[j], xi[j] + 1), length(uo_pars)), y = rep(uo_pars, each = 2), default.units = "native", id = rep(1:length(uo_pars), each = 2),
                          name = paste(lab, "item", j, "_uolines", sep = ""), gp = grid::gpar(col = uo_col, lwd = uo_lwd, lty = uo_lty))
          }
        }

        ## add box and axis
        grid::grid.rect(name = paste(lab, "_plot-box", sep = ""))
        grid::grid.xaxis(at = (xi[-(m+1)] + 0.5), label = names, main = TRUE, name = paste(lab, "_xaxis-bottom", sep = ""))
        grid::grid.yaxis(main = TRUE, name = paste(lab, "_yaxis-left", sep = ""))
        grid::upViewport()
        
        ## add descriptions
        grid::pushViewport(lmargin.vp)
        grid::grid.text(ylab, x = 0.2, rot = 90, name = paste(lab, "_ylab-left", sep = ""))
        grid::upViewport(2)
        
        ## go back to uper vp
        grid::upViewport(2)
    }

    ## return
    return(panelfun)
}
class(node_regionplot) <- "grapcon_generator"


## bradley-terry plot visualization function
node_btplot <- function(mobobj, id = TRUE, worth = TRUE, names = TRUE,
  abbreviate = TRUE, index = TRUE, ref = TRUE,col = "black", linecol = "lightgray",
  cex = 0.5, pch = 19, xscale = NULL, yscale = NULL, ylines = 1.5)
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
         
    ## panel function for bt plots in nodes
    rval <- function(node) {

      ## node index
      id <- id_node(node)
    
      ## dependent variable setup
      cfi <- cf[id,]

      ## viewport setup
      top_vp <- grid::viewport(layout = grid::grid.layout(nrow = 2, ncol = 3,
                                   widths = grid::unit(c(ylines, 1, 1), c("lines", "null", "lines")),  
                                   heights = grid::unit(c(1, 1), c("lines", "null"))),
                               width = grid::unit(1, "npc"), 
                               height = grid::unit(1, "npc") - grid::unit(2, "lines"),
                               name = paste("node_btplot", id, sep = ""))
      grid::pushViewport(top_vp)
      grid::grid.rect(gp = grid::gpar(fill = "white", col = 0))

      ## main title
      top <- grid::viewport(layout.pos.col = 2, layout.pos.row = 1)
      grid::pushViewport(top)
      mainlab <- paste(ifelse(id, paste("Node", id, "(n = "), ""),
        	       info_node(node)$nobs, ifelse(id, ")", ""), sep = "")
      grid::grid.text(mainlab)
      grid::popViewport()

      ## actual plot  
      plot_vpi <- grid::viewport(layout.pos.col = 2, layout.pos.row = 2,
        xscale = xscale, yscale = yscale, name = paste("node_btplot", id, "plot", sep = ""))
      grid::pushViewport(plot_vpi)

      grid::grid.lines(xscale, c(cf_ref, cf_ref), gp = grid::gpar(col = linecol), default.units = "native")
      if(index) {
        grid::grid.lines(x, cfi, gp = grid::gpar(col = col, lty = 2), default.units = "native")
        grid::grid.points(x, cfi, gp = grid::gpar(col = col, cex = cex), pch = pch, default.units = "native")
        grid::grid.xaxis(at = x, label = if(names) names(cfi) else x)
      } else {  	
        if(names) grid::grid.text(names(cfi), x = x, y = cfi, default.units = "native")
          else grid::grid.points(x, cfi, gp = grid::gpar(col = col, cex = cex), pch = pch, default.units = "native")
      }
      grid::grid.yaxis(at = c(ceiling(yscale[1] * 100)/100, floor(yscale[2] * 100)/100))
      grid::grid.rect(gp = grid::gpar(fill = "transparent"))

      grid::upViewport(2)
    }
	    
    return(rval)
}
class(node_btplot) <- "grapcon_generator"
