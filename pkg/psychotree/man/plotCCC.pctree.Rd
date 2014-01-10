\name{plotCCC.pctree}
\alias{plotCCC.pctree}


\title{Category Characteristic Curves for Rating Scale and Partial
  Credit Tree Models}


\description{

  Draws a matrix of category characteristic curves for a given number
  of nodes and items of a fitted rating scale or partial credit tree
  model.

}


\usage{

\S3method{plotCCC}{pctree}(object, nodes = NULL, items = NULL,
  byrow = TRUE, sum0 = TRUE, names = NULL, ylim = c(0, 1), xlim = NULL,
  col = NULL, n = 101, main = "Category Characteristic Curves",
  xlab = "Latent Trait", ylab = "Probability", \dots)

}


\arguments{
  
  \item{object}{an object of class \code{"mob"} based on rating scale
    models fitted by \code{\link{RSModel}} or partial credit models
    fitted by \code{\link{PCModel}}.}

  \item{nodes}{integer specifying the nodes for which the category
    characteristic curves of the items given in \code{items} should be
    drawn. If \code{NULL} (the default), category characteristic curves
    for all nodes are drawn.}

  \item{items}{numeric vector specifying the items for which category
    characteristic curves should be drawn. If \code{NULL} (the default),
    category characteristic curves for all items are drawn.}
  
  \item{byrow}{logical. If \code{TRUE} (the default) nodes are drawn in
    the rows and items in the columns. Otherwhise (\code{FALSE}) nodes
    are drawn in the columns and items in the rows.}

  \item{sum0}{logical. If \code{TRUE} (the default), the threshold
    parameters are normed to sum zero, otherwise (\code{FALSE}), the first
    threshold parameter is set to zero.}

  \item{names}{character vector of length \code{length(nodes)} x
    \code{length(items)} specifying the labels used in the stripes of the
    plot. If \code{NULL} (the default) these are generated.}

  \item{ylim}{y axis limits, typically zero and one.}

  \item{xlim}{x axis limits.}

  \item{col}{character vector specifying color names used for drawing the category
    characteristic curves.}

  \item{n}{numeric of length one specifying number of points to choose
    for plotting the category characteristic curves.}

  \item{main}{character, specifying the overall title of the plot.}

  \item{xlab}{label for the x axis.}

  \item{ylab}{label for the y axis.}

  \item{\dots}{further arguments passed to the various grid drawing
  functions used to generate the plot.}

}


\seealso{\code{\link{plotCCC}}, \code{\link{plotCCC.PCModel}}, \code{\link{node_ccc}}}


\examples{

## Verbal aggression data
data("VerbalAggression", package = "psychotools")

## fit partial credit tree model for the second other-to-blame situation
pct <- pctree(resp[, 7:12, drop = FALSE] ~ anger + gender, data =
VerbalAggression)

## default plot
plotCCC(pct)

## showing only item 1 and 4 for the two terminal nodes
plotCCC(pct, items = c(2, 4))

## showing the cccs of node 1 but with items in the rows
##  and specyfied x axis range.
plotCCC(pct, node = 1, byrow = FALSE, xlim = c(-5, 5))

}


\keyword{hplot}
