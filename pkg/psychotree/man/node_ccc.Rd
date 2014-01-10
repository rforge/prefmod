\name{node_ccc}
\alias{node_ccc}


\title{Panel-Generating Function for Visualizing Rating Scale and
  Partial Credit Tree Models}


\description{

  Panel-generating function for visualizing the threshold parameters from the
  nodes in rating scale and partial credit tree models.

}


\usage{

node_ccc(mobobj, names = NULL, sum0 = TRUE, ylim = NULL,
  ylab = "Latent trait", col = NULL, n = 101, ccc_lwd = 1.25)

}


\arguments{

  \item{mobobj}{an object of class \code{"mob"} based on rating scale models
    fitted by \code{\link{RSModel}} or partial credit models fitted by \code{\link{PCModel}}.}

  \item{names}{character vector containing x-axis labels for the
  category characteristic curves for each item. If set to \code{NULL}
  (the default), these are generated.}

  \item{sum0}{logical. If \code{TRUE} (the default), the threshold
    parameters are normed to sum zero, otherwise (\code{FALSE}), the
    first threshold parameter is set to zero.}

  \item{ylim}{y axis limits}

  \item{ylab}{label for the y axis.}

  \item{col}{character vector with color names for the various category
    characteristic curves of an item.}

  \item{n}{numeric of length one specifying the number of points to choose
    for plotting the category characteristic curves.}

  \item{ccc_lwd}{numeric of length one specifying the line width used
  for of the category characteristic curves}
  
}


\details{

  The panel-generating function \code{node_ccc} can be used via the
  argument \code{terminal_panel} for the \code{plot} method for
  \code{"rstree"} and \code{"pctree"} objects.

  In the former case (\code{"rstree"} objects), the parameters of the
  rating scale models in the terminal nodes are transformed to threshold
  parameters like in the partial credit model. Based on these, the
  category characteristic curves are drawn.

  Please note, there exists another high level function
  (\code{\link{plotCCC.pctree}}) which also allows plotting category
  characteristic curves. In contrast to \code{node_ccc}, the items and
  nodes for which category characteristic curves should be drawn, can be
  specified when using \code{\link{plotCCC.pctree}}, therefore allowing
  a more fine grained control over the output. In addition the
  characteristic curves of \code{\link{plotCCC.pctree}} are display
  horizontally, instead of the vertical arrangement here.

}


\value{

  A panel function which can be supplied to the \code{plot} method for
  \code{"mob"} objects.

}


\seealso{\code{\link{plotCCC.pctree}}, \code{\link{rstree}}, \code{\link{pctree}}}


\keyword{hplot}
