\name{plot.btreg}
\alias{plot.btreg}

\title{Visualization of Bradley-Terry Models}

\description{
  Flexible function for visualizing coefficients/worth parameters
  from a Bradley-Terry model.
}

\usage{
\method{plot}{btreg}(x, 
  worth = TRUE, index = TRUE, names = TRUE, ref = TRUE, abbreviate = FALSE,
  type = NULL, lty = NULL, xlab = "Objects", ylab = NULL, \dots)
}

\arguments{
  \item{x}{an object of class \code{"btreg"} as fitted by \code{\link{btreg}}.}
  \item{worth}{logical. Should worth parameters (or their logs) be visualized?}
  \item{index}{logical. Should different indexes for different stimuli be used?}
  \item{names}{logical. Should the names for the objects be displayed?}
  \item{ref}{logical. Should a horizontal line for the reference level be drawn?}
  \item{abbreviate}{logical or numeric. Should object names be abbreviated?
    If numeric this controls the length of the abbreviation.}
  \item{type, lty, xlab, ylab, \dots}{graphical parameters.}
}

\seealso{\code{\link{bttree}}}

\examples{
data("CEMSChoice", package = "prefmod2")
cems <- btreg(preference ~ 1, data = CEMSChoice)
plot(cems)
}

\keyword{hplot}
