\name{node_rsmplot}
\alias{node_rsmplot}


\title{Panel-Generating Function for Visualizing Rating Scale Tree Models}


\description{

  Panel-generating function for visualizing the item and threshold parameters
  from the nodes in rating scale tree models.

}


\usage{

node_rsmplot(mobobj, names.arg = NULL, ref = NULL,
  ylab = "Latent trait", ylim = NULL, col = c("gray", "gray"),
  pch = c(21, 22), lty = c(2, 3), refline = TRUE, reflinecol = "lightgray")

}


\arguments{

  \item{mobobj}{an object of class \code{"mob"} based on rating scale models
    fitted by \code{\link{RSModel}}.}

  \item{names.arg}{character vector containing x-axis labels for the item- and
    treshold parameters. If set to \code{NULL} (the default), these are generated.}

  \item{ref}{a vector of labels or position indices of item parameters
    which should be used as restriction/for normalization. If
    \code{NULL} (the default), all items are used (sum zero restriction).}

  \item{ylab}{label for the y axis.}

  \item{ylim}{y axis limits}

  \item{col}{character vector of length two with color names for filling
  the points of the item and threshold parameters.}

  \item{pch}{numeric vector of length two, specifying the symbols used
    for drawing the item and threshold parameters (see
    \code{\link[graphics]{points}} for possible values).}
  
  \item{lty}{numeric vector of length two, specifying the line type used
    for the connection lines between the various item and
    threshold parameters (see \code{\link[graphics]{lines}} for possible values).}
  
  \item{refline}{logical. Should a horizontal line for the reference level
    (mean item difficulty) be drawn?}

  \item{reflinecol}{line color for reference line (if \code{ref}).}
  
}


\details{

  The panel-generating function \code{node_rsmplot} is called by the \code{plot}
  method for \code{"rstree"} objects and does not have to be called by the
  user directly.

}


\value{

  A panel function which can be supplied to the \code{plot} method for
  \code{"mob"} objects.

}


\seealso{\code{\link{rstree}}}


\keyword{hplot}
