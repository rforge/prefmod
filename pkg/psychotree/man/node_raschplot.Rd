\name{node_raschplot}
\alias{node_raschplot}

\title{Panel-Generating Function for Visualizing Rasch Tree Models}

\description{
  Panel-generating function for visualizing the item parameters
  from the nodes in Rasch tree models.
}

\usage{
node_raschplot(mobobj, id = TRUE, difficulty = TRUE,
  center = TRUE, index = TRUE, names = NULL, abbreviate = FALSE, ref = TRUE,
  col = cbind("lightgray", "black"), refcol = "lightgray", linecol = "black",
  lty = 2, cex = 0.5, pch = cbind(19, 1), xscale = NULL, yscale = NULL,
  xaxis = TRUE, yaxis = TRUE, ylines = 1.5)
}

\arguments{
  \item{mobobj}{an object of class \code{"mob"} based on Rasch models
    fitted by \code{\link{RaschModel}}.}
  \item{id}{logical. Should the node ID be displayed?}
  \item{difficulty}{logical. Should item difficulty (or alternatively: easiness)
    parameters be displayed?}
  \item{center}{logical. Should the item parameters be centered?}
  \item{index}{logical. Should different indexes for different items be used?}
  \item{names}{logical. Should the names for the objects be displayed?}
  \item{abbreviate}{logical or numeric. Should object names be abbreviated?
    If numeric this controls the length of the abbreviation.}
  \item{ref}{logical. Should a horizontal line for the reference level be drawn?}
  \item{col, pch, cex}{graphical appearance of plotting symbols. Can be of the same
    length as number of items. Additionally \code{col} and \code{pch} can
    be matrices with two columns resulting in two symbols being overplotted.}
  \item{refcol}{line color for reference line (if \code{ref}).}
  \item{linecol}{line color.}
  \item{lty}{line type.}
  \item{xscale, yscale}{x and y axis limits.}
  \item{xaxis, yaxis}{logical. Should a axes be plotted?}
  \item{ylines}{numeric. Number of lines used for y-axis labels.}
}

\details{
  The panel-generating function \code{node_raschplot} is called by the \code{plot}
  method for \code{"raschtree"} objects and does not have to be called by the
  user directly.
}

\value{
  A panel function which can be supplied to the \code{plot} method for
  \code{"mob"} objects.
}

\seealso{\code{\link{raschtree}}}

\keyword{hplot}
