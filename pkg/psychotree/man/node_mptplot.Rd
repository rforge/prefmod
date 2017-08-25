\name{node_mptplot}
\alias{node_mptplot}

\title{Panel-Generating Function for Visualizing MPT Tree Models}

\description{
  Panel-generating function for visualizing the model parameters
  from the nodes in MPT tree models.
}

\usage{
node_mptplot(mobobj, id = TRUE,
  names = TRUE, abbreviate = TRUE, index = TRUE, ref = TRUE,
  col = "black", linecol = "lightgray", bg = "white", cex = 0.5, pch = 19,
  xscale = NULL, yscale = c(0, 1), ylines = 1.5)
}

\arguments{
  \item{mobobj}{an object of class \code{"mob"} based on MPT models
    fitted by \code{\link[psychotools]{mptmodel}}.}
  \item{id}{logical. Should the node ID be displayed?}
  \item{names}{logical or character. Should the names for the parameters be
    displayed? If character, this sets the names.}
  \item{abbreviate}{logical or numeric. Should parameter names be abbreviated?
    If numeric this controls the length of the abbreviation.}
  \item{index}{logical or character. Should different indexes for different
    parameters be used? If character, this controls the order of labels given
    in \code{names}.}
  \item{ref}{logical. Should a horizontal line for the reference level be
    drawn?}
  \item{col, cex, pch}{graphical appearance of plotting symbols.}
  \item{linecol}{line color for reference line (if \code{ref}).}
  \item{bg}{color for background filling.}
  \item{xscale, yscale}{x and y axis limits.}
  \item{ylines}{numeric. Number of lines used for y-axis labels.}
}

\details{
  The panel-generating function \code{node_mptplot} is called by the \code{plot}
  method for \code{"mpttree"} objects and does not have to be called by the
  user directly.
}

\value{
  A panel function which can be supplied to the \code{plot} method for
  \code{"mob"} objects.
}

\seealso{\code{\link{mpttree}}}

\keyword{hplot}
