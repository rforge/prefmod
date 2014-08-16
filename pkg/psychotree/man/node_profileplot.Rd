\name{node_profileplot}
\alias{node_profileplot}
\alias{node_raschplot}

\title{Panel-Generating Function for Visualizing Rasch,
  Rating Scale and Partial Credit Tree Models}

\description{
  Panel-generating function for visualizing the parameters
  from the nodes in Rasch, rating scale or partial credit tree models.
}

\usage{
node_profileplot(mobobj, what = c("items", "thresholds", "discriminations"),
  paramarg = list(type = NULL, ref = NULL, alias = TRUE), id = TRUE,
  names = TRUE, abbreviate = TRUE, index = TRUE, ref = TRUE, col = "black",
  linecol = "lightgray", cex = 0.5, pch = 19, xscale = NULL, yscale =
  NULL, ylines = 1.5, ...)
}

\arguments{
  \item{mobobj}{an object of class \code{"mob"} based on rasch models
    fitted by \code{\link[psychotools]{raschmodel}}, rating scale models
    fitted by \code{\link[psychotools]{rsmodel}} or partial credit
    models fitted by \code{\link[psychotools]{pcmodel}}.}
  \item{what}{character, specifying the type of parameters to be plotted.}
  \item{paramarg}{list of arguments passed over to internal calls of
    \code{\link[psychotools]{itempar}}, \code{\link[psychotools]{threshpar}}
    or \code{\link[psychotools]{discrpar}}.}
  \item{id}{logical. Should the node ID be displayed?}
  \item{names}{logical or character. If \code{TRUE}, the names of the
    items are displayed. If a vector of characters, these are used as labels}
  \item{abbreviate}{logical. Should object names be abbreviated?
    If numeric this controls the length of the abbreviation.}
  \item{index}{logical. Should different indexes for different items be used?}
  \item{ref}{logical. Should a horizontal line for the reference level be drawn?}
  \item{col, pch, cex}{graphical appearance of plotting symbols.}
  \item{linecol}{character, specifying the line color to use for the
    reference line.}
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
