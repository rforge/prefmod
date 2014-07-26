### high-level convenience interface to mob()
pctree <- function (formula, data, na.action = na.pass, nullcats = c("keep", "downcode", "ignore"),
  reltol = 1e-10, deriv = c("sum", "diff"), maxit = 100L, ...)
{
  ## keep call
  cl <- match.call(expand.dots = TRUE)

  ## use dots for setting up mob_control
  control <- mob_control(...)
  control$ytype <- "matrix"

  ## control options for pcmfit
  pcmcontrol <- list(nullcats = nullcats, reltol = reltol, deriv = deriv, maxit = maxit)

  ## call mob
  m <- match.call(expand.dots = FALSE)
  m$fit <- pcmfit
  m$control <- control
  for(n in names(pcmcontrol)) if(!is.null(pcmcontrol[[n]])) m[[n]] <- pcmcontrol[[n]]
  if("..." %in% names(m)) m[["..."]] <- NULL
  m[[1L]] <- as.name("mob")
  rval <- eval(m, parent.frame())
  
  ## extend class and keep original call
  rval$info$call <- cl
  class(rval) <- c("pctree", class(rval))
  return(rval)
}

## glue code for calling pcmodel()
pcmfit <- function(y, x = NULL, start = NULL, weights = NULL, offset = NULL, ...,
  estfun = FALSE, object = FALSE)
{
  if(!(is.null(x) || NCOL(x) == 0L)) warning("x not used")
  if(!is.null(offset)) warning("offset not used")
  rval <- pcmodel(y, weights = weights, start = start, ..., hessian = object | estfun)
  rval <- list(
    coefficients = rval$coefficients,
    objfun = -rval$loglik,
    estfun = if(estfun) estfun.pcmodel(rval) else NULL,
    object = if(object) rval else NULL
  )
  return(rval)
}

## methods
print.pctree <- function(x,
  title = "Partial credit tree", objfun = "negative log-likelihood", ...)
{
  partykit::print.modelparty(x, title = title, objfun = objfun, ...)
}

predict.pctree <- function(object, newdata = NULL,
  type = c("node"), ...)
{
  ## type of prediction
  type <- match.arg(type)
  
  ## nodes can be handled directly
  if(type == "node") return(partykit::predict.modelparty(object, newdata = newdata, type = "node", ...))
  
  ## get default newdata otherwise
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

## worth.pctree <- function(object, node = NULL, ...)
## {
##   ids <- if(is.null(node)) nodeids(object, terminal = TRUE) else node
##   if(length(ids) == 1L) {
##     apply_to_models(object, node = ids, FUN = worth, drop = TRUE)
##   } else {
##     do.call("rbind", apply_to_models(object, node = ids, FUN = worth, drop = FALSE))
##   } 
## }

plot.pctree <- function(x, terminal_panel = node_effects,
  tp_args = list(), tnex = NULL, drop_terminal = NULL, ...)
{
  if(is.null(tnex)) tnex <- if(is.null(terminal_panel)) 1L else 2L
  if(is.null(drop_terminal)) drop_terminal <- !is.null(terminal_panel)
  partykit::plot.modelparty(x, terminal_panel = terminal_panel,
    tp_args = tp_args, tnex = tnex, drop_terminal = drop_terminal, ...)
}

## visualization function
node_effects <- function(mobobj, names = NULL, type = c("mode", "median", "mean"),
                         ref = NULL, ylab = "Latent trait", ylim = NULL, off = 0.1, col_fun = gray.colors,
                         uo_show = TRUE, uo_col = "red", uo_lty = 2, uo_lwd = 1.25)    
{
    ## check input
    stopifnot(!is.null(mobobj))
    stopifnot(off >= 0)
    type <- match.arg(type)

    ## function to extract absolute item threshold parameters from model objects in terminal nodes
    threshparlst <- function (node) coef(threshpar(node, ref = ref, type = type, alias = TRUE, relative = FALSE, cumulative = FALSE, vcov = FALSE), type = "list")

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
        top.vp <- viewport(layout = grid.layout(nrow = 2, ncol = 1, widths = unit(1, "null"), heights = unit(c(1, 1), c("lines", "null"))),
                           width = unit(1, "npc"), height = unit(1, "npc") - unit(2, "lines"), name = paste(lab, "_effects", sep = ""))
        pushViewport(top.vp)
        grid.rect(gp = gpar(fill = "white", col = 0), name = paste(lab, "_border", sep = ""))

        ## main title
        pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 1, name = paste(lab, "_title_vp", sep = "")))
        grid.text(paste("Node ", id, " (n = ", info_node(node)$nobs, ")", sep = ""), name = paste(lab, "_title", sep = ""))
        upViewport()
        
        ## finally the actual plot
        pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 2, name = lab))
        lab <- paste(lab, "_plot", sep = "")

        ## setup plotting area (3x3)
        wcol <- if (is.null(ylab)) c(2.5, 1, 1) else c(4, 1, 1)
        hrow <- c(0.5, 1, 1)

        top.vp <- viewport(layout=grid.layout(nrow = 3, ncol = 3, widths = unit(wcol, c("lines", "null", "lines")), heights = unit(hrow, c("lines", "null", "lines"))), name = paste(lab, "_top_vp", sep = ""))
        bmargin.vp <- viewport(layout.pos.row = 3, layout.pos.col = 2, name = paste(lab, "_bottom-margin_vp", sep = ""))
        lmargin.vp <- viewport(layout.pos.row = 2, layout.pos.col = 1, name = paste(lab, "_left-margin_vp", sep = ""))
        rmargin.vp <- viewport(layout.pos.row = 2, layout.pos.col = 3, name = paste(lab, "_right-margin_vp", sep = ""))
        plot.vp <- viewport(layout.pos.row = 2, layout.pos.col = 2, name = paste(lab, "_vp", sep = ""), xscale = xlim, yscale = ylim)
        pushViewport(top.vp)
        pushViewport(plot.vp)

        ## plot rectangles per item
        for (j in seq_along(delta_sorted)) {
            ncat <- length(delta_sorted[[j]]) + 1
            grid.rect(x = rep.int(xi[j], ncat), y = c(ylim[1], delta_sorted[[j]]), width = rep.int(1, ncat),
                      height = diff.default(c(ylim[1], delta_sorted[[j]], ylim[2])), just = c("left", "bottom"),
                      gp = gpar(fill = col_fun(ncat)), default.units = "native", name = paste(lab, "_item", j, "_rect", sep = ""))
        }

        ## if requested: indicate unordered parameters
        if (uo_show && type == "mode") {
          uo_items <- which(!sapply(mapply(all.equal, delta_sorted, delta_unsorted, check.attributes = FALSE, SIMPLIFY = FALSE, USE.NAMES = FALSE), is.logical))
          for (j in uo_items) {
            uo_pars <- setdiff(delta_unsorted[[j]], delta_sorted[[j]])
            grid.polyline(x = rep(c(xi[j], xi[j] + 1), length(uo_pars)), y = rep(uo_pars, each = 2), default.units = "native", id = rep(1:length(uo_pars), each = 2),
                          name = paste(lab, "item", j, "_uolines", sep = ""), gp = gpar(col = uo_col, lwd = uo_lwd, lty = uo_lty))
          }
        }

        ## add box and axis
        grid.rect(name = paste(lab, "_plot-box", sep = ""))
        grid.xaxis(at = (xi[-(m+1)] + 0.5), label = names, main = TRUE, name = paste(lab, "_xaxis-bottom", sep = ""))
        grid.yaxis(main = TRUE, name = paste(lab, "_yaxis-left", sep = ""))
        upViewport()
        
        ## add descriptions
        pushViewport(lmargin.vp)
        grid.text(ylab, x = 0.2, rot = 90, name = paste(lab, "_ylab-left", sep = ""))
        upViewport(2)
        
        ## go back to uper vp
        upViewport(2)
    }

    ## return
    return(panelfun)
}
class(node_effects) <- "grapcon_generator"
