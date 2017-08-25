\name{node_profileplot}
\alias{node_profileplot}
\alias{node_raschplot}

\title{Panel-Generating Function for Visualizing IRT Tree Models}

\description{
  Panel-generating function for visualizing the parameters
  from the nodes in Rasch, rating scale or partial credit tree models.
}

\usage{
node_profileplot(mobobj, what = c("items", "thresholds", "discriminations"),
  parg = list(type = NULL, ref = NULL, alias = TRUE), id = TRUE,
  names = FALSE, abbreviate = TRUE, index = TRUE, ref = TRUE, col = "black",
  border = col, linecol = "black", refcol = "lightgray", bg = "white",
  cex = 0.5, pch = 21, xscale = NULL, yscale = NULL, ylines = 2, ...)
}

\arguments{
  \item{mobobj}{an object of class \code{"mob"} based on rasch models
    fitted by \code{\link[psychotools]{raschmodel}}, rating scale models
    fitted by \code{\link[psychotools]{rsmodel}} or partial credit
    models fitted by \code{\link[psychotools]{pcmodel}}.}
  \item{what}{character, specifying the type of parameters to be plotted.}
  \item{parg}{list of arguments passed over to internal calls of
    \code{\link[psychotools]{itempar}}, \code{\link[psychotools]{threshpar}}
    or \code{\link[psychotools]{discrpar}}.}
  \item{id}{logical. Should the node ID be displayed?}
  \item{names}{logical or character. If \code{TRUE}, the names of 
    the items are displayed on the x-axis. If \code{FALSE}, numbers of
    items are shown. Alternatively a character vector of the same
    length as the number of items can be supplied.}
  \item{abbreviate}{logical. Should item names be abbreviated?
    If numeric this controls the length of the abbreviation.}
  \item{index}{logical. Should different indexes for different items be used?}
  \item{ref}{logical. Should a horizontal line for the reference level be drawn?}
  \item{col, border, pch, cex}{graphical appearance of plotting symbols.}
  \item{linecol, refcol}{character, specifying the line color to use for the
    profile lines and reference line, respectively.}
  \item{bg}{color for background filling.}
  \item{xscale, yscale}{x and y axis limits.}
  \item{ylines}{numeric. Number of lines used for y-axis labels.}
  \item{\dots}{further arguments currently not used.}
}

\details{
  The panel-generating function \code{node_profileplot} is called by the
  \code{plot} method for \code{"raschtree"} objects and does not have to
  be called by the user directly.

  Prior to version 0.13-1, \code{node_profileplot}
  was called \code{node_raschplot}.
}

\value{
  A panel function which can be supplied to the \code{plot} method for
  \code{"mob"} objects.
}

\seealso{\code{\link{raschtree}}}

\keyword{hplot}
