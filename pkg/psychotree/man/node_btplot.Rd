\name{node_btplot}
\alias{node_btplot}

\title{Panel-Generating Function for Visualizing Bradley-Terry Tree Models}

\description{
  Panel-generating function for visualizing the worth parameters
  from the nodes in Bradley-Terry tree models.
}

\usage{
node_btplot(mobobj, id = TRUE,
  worth = TRUE, names = TRUE, abbreviate = TRUE, index = TRUE, ref = TRUE,
  col = "black", linecol = "lightgray", cex = 0.5, pch = 19,
  xscale = NULL, yscale = NULL, ylines = 1.5)
}

\arguments{
  \item{mobobj}{an object of class \code{"mob"} based on Bradley-Terry models
    fitted by \code{\link[psychotools]{btmodel}}.}
  \item{id}{logical. Should the node ID be displayed?}
  \item{worth}{logical. Should worth parameters (or their logs) be visualized?}
  \item{names}{logical. Should the names for the objects be displayed?}
  \item{abbreviate}{logical or numeric. Should object names be abbreviated?
    If numeric this controls the length of the abbreviation.}
  \item{index}{logical. Should different indexes for different stimuli be used?}
  \item{ref}{logical. Should a horizontal line for the reference level be drawn?}
  \item{col, cex, pch}{graphical appearance of plotting symbols.}
  \item{linecol}{line color for reference line (if \code{ref}).}
  \item{xscale, yscale}{x and y axis limits.}
  \item{ylines}{numeric. Number of lines used for y-axis labels.}
}

\details{
  The panel-generating function \code{node_btplot} is called by the \code{plot}
  method for \code{"bttree"} objects and does not have to be called by the
  user directly.
}

\value{
  A panel function which can be supplied to the \code{plot} method for
  \code{"mob"} objects.
}

\seealso{\code{\link{bttree}}}

\keyword{hplot}
