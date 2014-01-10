\name{node_effects}
\alias{node_effects}


\title{Panel-Generating Function for Visualizing Rating Scale and
  Partial Credit Tree Models}


\description{

  Panel-generating function for visualizing the threshold parameters from the
  nodes in rating scale and partial credit tree models.

}


\usage{

node_effects(mobobj, names = NULL, type = c("mode", "median", "mean"),
  ref = NULL, ylab = "Latent trait", ylim = NULL, off = 0.1, col_fun = gray.colors,
  uo_show = TRUE, uo_col = "red", uo_lty = 2, uo_lwd = 1.25, ...)

}


\arguments{

  \item{mobobj}{an object of class \code{"mob"} based on rating scale models
    fitted by \code{\link{RSModel}} or partial credit models fitted by \code{\link{PCModel}}.}

  \item{names}{character vector containing x-axis labels for the
    threshold parameters. If set to \code{NULL} (the default), these are
    generated.}
  
  \item{type}{character, specifying which type of threshold parameters
    are to be used to mark the category regions per item in the plot
    (see \code{\link[psychotools]{plot.PCModel}} for details).}
  
  \item{ref}{a vector of labels or position indices of item parameters
    which should be used as restriction/for normalization. If
    \code{NULL} (the default), all items are used (sum zero restriction).}

  \item{ylim}{y axis limits}

  \item{ylab}{label for the y axis.}

  \item{off}{the distance (in scale units) between two item rectangles.}

  \item{col_fun}{function. Function to use for creating the color
  palettes for the rectangles. Per default \code{gray.colors} is
  used. Be aware that \code{col_fun} should accept as first argument an
  integer specifying the number of colors to create.}

  \item{uo_show}{logical. If set to \code{TRUE} (the default),
    disordered threshold parameters are indicated by a horizontal line
    (if \code{cut} is set to \code{"modus"}).}

  \item{uo_col}{character, color of indication lines (if \code{uo_show}).}

  \item{uo_lty}{numeric, line typ of indication lines (if
    \code{uo_show}).}

  \item{uo_lwd}{numeric, line width of indication lines (if
    \code{uo_show}).}

  \item{...}{additional arguments passed to col_fun.}

}


\details{

  The panel-generating function \code{node_effects} is called by the
  \code{plot} method for \code{"rstree"} and \code{"pctree"} objects by
  default and does not have to be called by the user directly. 

  See \code{\link[psychotools]{plot.PCModel}} for details and references
  of the drawn effect plots and possible values and their meaning for
  argument \code{cut} (taken by \code{node_effects}).

}


\value{

  A panel function which can be supplied to the \code{plot} method for
  \code{"mob"} objects.

}


\seealso{\code{\link{rstree}}, \code{\link{pctree}},
\code{\link[psychotools]{plot.PCModel}}}


\keyword{hplot}
