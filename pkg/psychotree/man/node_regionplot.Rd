\name{node_regionplot}
\alias{node_regionplot}
\alias{node_effects}

\title{Panel-Generating Function for Visualizing Rasch, Rating Scale and
  Partial Credit Tree Models}

\description{
  Panel-generating function for visualizing the absolute item threshold
  parameters from the nodes in rasch, rating scale and partial credit
  tree models.
}

\usage{
node_regionplot(mobobj, names = FALSE, abbreviate = TRUE,
  type = c("mode", "median", "mean"), ref = NULL, ylim = NULL,
  off = 0.1, col_fun = gray.colors, uo_show = TRUE, uo_col = "red",
  uo_lty = 2, uo_lwd = 1.25, ylines = 2)
}

\arguments{
  \item{mobobj}{an object of class \code{"mob"} based on rasch models
    fitted by \code{\link[psychotools]{raschmodel}}, rating scale models
    fitted by \code{\link[psychotools]{rsmodel}} or partial credit
    models fitted by \code{\link[psychotools]{pcmodel}}.}
  \item{names}{logical or character. If \code{TRUE}, the names of 
    the items are displayed on the x-axis. If \code{FALSE}, numbers of
    items are shown. Alternatively a character vector of the same
    length as the number of items can be supplied.}
  \item{abbreviate}{logical. Should item names be abbreviated?
    If numeric this controls the length of the abbreviation.}
  \item{type}{character, specifying which type of threshold parameters
    are to be used to mark the category regions per item in the plot
    (see \code{\link[psychotools]{regionplot}} for details).}
  \item{ref}{a vector of labels or position indices of item parameters
    which should be used as restriction/for normalization. If
    \code{NULL} (the default), all items are used (sum zero
    restriction). See \code{\link[psychotools]{threshpar}} for more details.}
  \item{ylim}{y axis limits}
  \item{off}{numeric, the distance (in scale units) between two item rectangles.}
  \item{col_fun}{function. Function to use for creating the color
    palettes for the rectangles. Per default \code{gray.colors} is
    used. Be aware that \code{col_fun} should accept as first argument an
    integer specifying the number of colors to create.}
  \item{uo_show}{logical. If set to \code{TRUE} (the default),
    disordered absolute item threshold parameters are indicated by a
    horizontal line (only if \code{type} is set to \code{"mode"}).}
  \item{uo_col}{character, color of indication lines (if \code{uo_show}).}
  \item{uo_lty}{numeric. Line typ of indication lines (if \code{uo_show}).}
  \item{uo_lwd}{numeric. Line width of indication lines (if \code{uo_show}).}
  \item{ylines}{numeric. Number of lines used for y-axis labels.}
}

\details{
  The panel-generating function \code{node_regionplot} is called by the
  \code{plot} method of \code{"rstree"} and \code{"pctree"} objects by
  default and does not have to be called by the user directly.

  See \code{\link[psychotools]{regionplot}} for details and references
  of the drawn region plots and possible values and their meaning for
  the argument \code{type} (taken by \code{node_regionplot}).

  Prior to version 0.13-1, \code{node_regionplot}
  was called \code{node_effects}.
}

\value{
  A panel function which can be supplied to the \code{plot} method
  for \code{"mob"} objects.
}

\seealso{\code{\link{raschtree}}, \code{\link{rstree}},
  \code{\link{pctree}}, \code{\link[psychotools]{regionplot}}}

\keyword{hplot}
