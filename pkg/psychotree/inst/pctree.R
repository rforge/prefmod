if(FALSE) {
library("partykit")
library("psychotools")
source("../R/pctree-psychotools.R")
data("DIFSim", package = "psychotree")
mb <-  mob(resp ~ age + gender + motivation, data = DIFSim, fit = pcmfit, control = mob_control(ytype = "matrix"))
pct <- pctree(resp ~ age + gender + motivation, data = DIFSim)
}

### high-level convenience interface to mob()
pctree <- function (formula, data, na.action = na.pass, nullcats = c("keep", "downcode", "ignore"),
                    reltol = 1e-10, deriv = c("sum", "diff"), ...)
{
  ## keep call
  cl <- match.call(expand.dots = TRUE)

  ## use dots for setting up mob_control
  control <- mob_control(...)
  control$ytype <- "matrix"

  ## control options for pcmfit
  pcmcontrol <- list(nullcats = nullcats, reltol = reltol, deriv = deriv)

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

## glue code for calling PCModel.fit()
pcmfit <- function(y, x = NULL, start = NULL, weights = NULL, offset = NULL, ...,
  estfun = FALSE, object = FALSE)
{
  if(!(is.null(x) || NCOL(x) == 0L)) warning("x not used")
  if(!is.null(offset)) warning("offset not used")
  rval <- PCModel.fit(y, weights = weights, start = start, ..., hessian = object | estfun)
  rval <- list(
    coefficients = rval$coefficients,
    objfun = -rval$loglik,
    estfun = if(estfun) estfun.PCModel(rval) else NULL,
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

## rewrite of apply_to_models() -- "..." is now passed to FUN
apply_to_models2 <- function(object, node = NULL, FUN = NULL, drop = FALSE, ...) {
  if(is.null(node)) node <- nodeids(object, terminal = FALSE)
  if(is.null(FUN)) FUN <- function(object) object  
  rval <- if("object" %in% object$info$control$terminal) {
    nodeapply(object, node, function(n) FUN(info_node(n)$object, ...))
  } else {
    lapply(refit.modelparty(object, node, drop = FALSE), FUN, ...)
  }
  names(rval) <- node
  if(drop & length(node) == 1L) rval <- rval[[1L]]
  return(rval)
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

    ## FIXME: glue code as long as threshold() and itempar() is not exported:
    ## get one full model, determine class, set appropriate threshold and itempar function
    model <- apply_to_models2(mobobj, node = 1L, FUN = NULL)[[1]]
    RSM <- FALSE
    if (inherits(model[[1]], "RaschModel")) {
        threshold <- threshold.RaschModel
        itempar <- function (node) lapply(itempar.RaschModel(node, ref = ref, vcov = FALSE, simplify = FALSE), function (j) diff.default(c(0, j)))
    } else if (inherits(model[[1]], "RSModel")) {
        threshold <- threshold.RSModel
        itempar <- function (node) lapply(as.list(itempar.RSModel(node, ref = ref, vcov = FALSE, simplify = FALSE)[[1]]),
                                          function (beta) diff(0:length(ip[[2]]) * beta + c(0, ip[[2]])))
    } else {
        threshold <- threshold.PCModel
        itempar <- function (node) lapply(itempar.PCModel(node, ref = ref, vcov = FALSE, simplify = FALSE), function (j) diff.default(c(0, j)))
    }

    ## setup threshold parameters
    node <- nodeids(mobobj, terminal = TRUE)
    delta_lst <- apply_to_models2(mobobj, node, FUN = threshold, type = type, ref = ref, simplify = FALSE)
    
    ## if requested and type = 'mode' check for unordered thresholds
    if (uo_show && type == "mode") {
        ip_lst <- apply_to_models2(mobobj, node, FUN = itempar)
        names(ip_lst) <- node
    }

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
        cf <- delta_lst[[id]]
        lab <- paste("node", id, sep = "")
        
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
        for (j in seq_along(cf)) {
            
            ncat <- length(cf[[j]]) + 1
            grid.rect(x = rep.int(xi[j], ncat), y = c(ylim[1], cf[[j]]), width = rep.int(1, ncat), height = diff.default(c(ylim[1], cf[[j]], ylim[2])),
                      just = c("left", "bottom"), gp = gpar(fill = col_fun(ncat)), default.units = "native",
                      name = paste(lab, "_item", j, "_rect", sep = ""))

        }

        ## if requested: indicate unordered parameters
        if (uo_show && type == "mode") {
            ip <- ip_lst[[id]]
            uo_items <- which(!sapply(mapply(all.equal, cf, ip, check.attributes = FALSE, SIMPLIFY = FALSE, USE.NAMES = FALSE), is.logical))
            for (j in uo_items) {
                uo_pars <- setdiff(ip[[j]], cf[[j]])
                grid.polyline(x = rep(c(xi[j], xi[j] + 1), length(uo_pars)), y = rep(uo_pars, each = 2), 
                              default.units = "native", id = rep(1:length(uo_pars), each = 2),
                              name = paste(lab, "item", j, "_uolines", sep = ""),
                              gp = gpar(col = uo_col, lwd = uo_lwd, lty = uo_lty))
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

bread.PCModel <- function(x, ...) x$vcov * x$n

estfun.PCModel <- function (x, ...) {
  ## get relevant informations
  dat <- x$data                    # completely cleaned (downcoded, null cats treatment, weights) data.
  weights_org <- weights(x)
  weights <- weights_org[weights_org > 0]
  n <- nrow(dat)
  m <- ncol(dat)
  oj_vec <- x$categories
  oj <- sapply(x$categories, length)
  npar_all <- sum(oj)
  npar_ident <- x$df
  ptot <- rowSums(dat, na.rm = TRUE) + 1 # +1 because gamma of score 0 is in row 1.

  ## helper variables
  parindex <- unlist(oj_vec)
  itemindex <- rep.int(1:m, oj)

  ## calculate gradient
  if (!x$na) {

    ## select gamma zero and first derivatives with ptot
    gamma0 <- x$esf[[1]][ptot]
    gamma1 <- apply(x$esf[[2]], 2, "[", ptot)

    ## construct data matrix ('selection' matrix, 0/1, cols = parameters)
    if (!is.null(x$nullcats)) { ## null cats & strategy == 'keep', remove column of unidentified par
      est_par <- !unlist(x$nullcats)
      gamma1 <- gamma1[, est_par, drop = FALSE]
    }
    xmat <- matrix(FALSE, nrow = n, ncol = npar_all)
    for (i in 1:n) xmat[i, ] <- dat[i, itemindex] == parindex 

    ## calculate gradient
    agrad <- weights * (- xmat + (gamma1 / gamma0))

  } else {

    ## return value & helper variables
    agrad <- matrix(0, nrow = n, ncol = npar_all)
    mv <- 1:m

    ## observed NA patterns 
    na_patterns <- factor(apply(is.na(dat), 1, function(z) paste(which(z), collapse = "\r")))

    ## loop through na patterns, select derivatives and calculate gradient
    for(i in seq_len(nlevels(na_patterns))) {

      ## parse NA patterns and setup necessary stuff for gradient calculation
      lev_i <- levels(na_patterns)[i]
      na_i <- which(na_patterns == lev_i)
      n_na_i <- length(na_i)
      mv_i <- as.integer(strsplit(lev_i, "\r")[[1]])
      mv_i <- if(length(mv_i) < 1) mv else mv[-mv_i]
      oj_i <- oj[mv_i]
      oj_vec_i <- oj_vec[mv_i]
      weights_i <- weights[na_i]
      ptot_i <- ptot[na_i]
      dat_i <- dat[na_i, , drop = FALSE]

      ## select gamma zero and first derivatives with ptot_i
      gamma0_i <- x$esf[[i]][[1]][ptot_i]
      gamma1_i <- apply(x$esf[[i]][[2]], 2, "[", ptot_i)
      if (!is.matrix(gamma1_i)) gamma1_i <- matrix(gamma1_i, nrow = 1)

      ## construct data matrix ('selection' matrix, 0/1, cols = parameters) for NA group i
      parindex_i <- unlist(oj_vec_i)
      itemindex_i <- rep.int(mv_i, oj_i)

      if (!is.null(x$nullcats)) { ## null cats & strategy == 'keep', remove column of unidentified par
        est_par_i<- !unlist(x$nullcat[mv_i])
        gamma1_i <- gamma1_i[, est_par_i, drop = FALSE]
      }

      xmat_i <- matrix(FALSE, nrow = n_na_i, ncol = sum(oj_i))
      for (i in 1:n_na_i) xmat_i[i, ] <- dat_i[i, itemindex_i] == parindex_i

      ## finally: the gradient for NA group i
      agrad[na_i, itemindex %in% mv_i] <- weights_i * (- xmat_i + (gamma1_i/ gamma0_i))
    }

  }  

  ## collect and return matrix of initial size with gradients plugged in.
  grad <- matrix(0, ncol = npar_ident, nrow = length(weights_org))
  grad[weights_org > 0, ] <- agrad[, -1, drop = FALSE]
  return(grad)
}
