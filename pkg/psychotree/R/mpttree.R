# TODO:
# - mpttree()
#   * formula interface only here! DONE
#   * a, b, c extraction DONE
#   * treeid processing: necessary?


## High-level convenience interface
mpttree <- function(formula, data, treeid = NULL, na.action = na.pass,
                    minsplit = 5, ...)
{
  ## Transform formula
  stopifnot(length(formula) > 2)
  formula <- formula(terms(formula, data = data))
  ff <- y ~ 1 | x
  ff[[2]] <- formula[[2]]
  ff[[3]][[3]] <- formula[[3]][[3]]  # keep partitioning vars only
  mptform <- formula[[3]][[2]]       # model structure
  tnames <- all.vars(mptform)        # names(theta)

  ## Extract model structure from formula
  aa <- bb <- array(NA,
    c(max(sapply(gregexpr("\\+", mptform), function(x) sum(x > 0))) + 1,
      length(mptform) - 1, length(tnames)))
  cc <- matrix(1, dim(aa)[1], dim(aa)[2])

  terms <- strsplit(as.character(mptform[-1]), "\\+")
  terms <- lapply(terms, function(x) gsub(" ", "", x))    # remove white space

  for(j in 1:dim(aa)[2]){
    for(i in 1:sapply(terms, length)[j]){
      pterms <- strsplit(terms[[j]][i], "\\*")[[1]]
      cval <- sum(as.numeric(grep("^[0-9]+$", pterms, value=TRUE)))
      if(cval > 0) cc[i,j] <- cval
      for(s in seq_along(tnames)){
        tname <- tnames[s]

        aa[i,j,s] <- sum(grepl(paste0("^", tname, "$"), pterms))
        powix <- grepl(paste0("^", tname, "\\^[0-9]+"), pterms)
        aa[i,j,s] <- sum(aa[i,j,s],
          as.numeric(gsub(paste0("^", tname, "\\^([0-9]+)"), "\\1",
            pterms)[powix]))

        ## Brackets () are optional
        bb[i,j,s] <- sum(grepl(paste0("^\\(?1-", tname, "\\)?$"), pterms))
        powix <- grepl(paste0("^\\(1-", tname, "\\)\\^[0-9]+"), pterms)
        bb[i,j,s] <- sum(bb[i,j,s],
          as.numeric(gsub(paste0("^\\(1-", tname, "\\)\\^([0-9]+)"), "\\1",
            pterms)[powix]))
      }
    }
  }
  dimnames(aa)[[3]] <- dimnames(bb)[[3]] <- as.list(tnames)

  ## call mob()
  rval <- mob(ff, data = data,
              model = mptModel(treeid = treeid, mptform = mptform,
                               mptstruc = list(a=aa, b=bb, c=cc)),
    control = mob_control(minsplit = minsplit, ...), na.action = na.action)

  ## add class and return
  structure(list(mob = rval), class = "mpttree")
}


## Convenience plotting
plot.mpttree <- function(x, terminal_panel = node_mptplot, tnex = 2, ...) {
  plot(x$mob, terminal_panel = terminal_panel, tnex = tnex, tp_args = list(...))
}


## Hand-crafted "Next()" to bridge to
## un-exported S4 classes "mob"/"BinaryTree"
deviance.mpttree <- function(object, ...) deviance(object$mob, ...)
logLik.mpttree <- function(object, ...) logLik(object$mob, ...)
sctest.mpttree <- function(x, ...) sctest(x$mob, ...)
weights.mpttree <- function(object, ...) weights(object$mob, ...)
summary.mpttree <- function(object, ...) summary(object$mob, ...)
print.mpttree <- function(x, ...) {
  print(x$mob, ...)
  invisible(x)
}


## Parameters for MPT trees
coef.mpttree <- function (object, node = NULL, ...)
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


## Visualization function
node_mptplot <- function(mobobj, id = TRUE,
  names = TRUE, abbreviate = TRUE, index = TRUE, ref = TRUE,
  col = "black", linecol = "lightgray", cex = 0.5, pch = 19, xscale = NULL,
  yscale = c(0, 1), ylines = 1.5)
{
    ## extract parameter of interest
    node <- 1:max(party:::terminal_nodeIDs(mobobj@tree))
    cf <- t(sapply(nodes(mobobj, node),
      function(z) coef(z$model, all = FALSE, ref = TRUE)))
    rownames(cf) <- node

    ## reference
    if(ref)
      cf_ref <- 1/2

    ## labeling
    if(is.character(names)) {
      colnames(cf) <- names
      names <- TRUE
    }
    if(is.character(index)) {
      cf <- cf[, index]  # reorder labels
      index <- TRUE
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

    ## panel function for mpt plots in nodes
    rval <- function(node) {

        ## dependent variable setup
        cfi <- cf[node$nodeID,]

        ## viewport setup
        top_vp <- viewport(layout = grid.layout(nrow = 2, ncol = 3,
                           widths = unit(c(ylines, 1, 1), c("lines", "null", "lines")),
                           heights = unit(c(1, 1), c("lines", "null"))),
                           width = unit(1, "npc"),
                           height = unit(1, "npc") - unit(2, "lines"),
                           name = paste("node_mptplot", node$nodeID, sep = ""))
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
            name = paste("node_mptplot", node$nodeID, "plot", sep = ""))
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
class(node_mptplot) <- "grapcon_generator"

