\name{plot.paircomp}
\alias{plot.paircomp}

\title{Plotting Paired Comparison Data}

\description{
  Plotting the frequency table from \code{"paircomp"} data.
}

\usage{
  \method{plot}{paircomp}(x, off = 0.05,
  xlab = "Proportion of comparisons", ylab = "", tol.xlab = 0.05,
  abbreviate = TRUE, hue = NULL, chroma = 40, luminance = 80)
}

\arguments{
  \item{x}{an object of class \code{"paircomp"}.}
  \item{off}{numeric. Offset between segments on the y-axis.}
  \item{xlab, ylab}{character. Axis labels.}
  \item{tol.xlab}{numeric. convenience tolerance parameter for x-axis annotation.
    If the distance between two labels drops under this threshold, they
    are plotted equidistantly.}
  \item{abbreviate}{logical or integer. Should object labels be abbreviated?
    Alternative an integer with the desired abbreviation length. The default
    is some heuristic based on the length of the labels.}
  \item{hue}{numeric. A vector of hues in [0, 360], recycled to the
    number of objects compared in \code{x}. A sequential palette
    is computed for each hue, see below.}
  \item{chroma}{numeric. Maximum chroma in the palette.}
  \item{luminance}{numeric. Minimum (and maximum) luminance in the
    palette. If omitted, the maximum is set to 95.}
}

\details{
  The \code{plot} method creates a frequency table (using \code{summary})
  and visualizes this using a sort of spine plot with HCL-based
  diverging palettes. See Zeileis, Hornik, Murrell (2009) for the
  underlying ideas.
}

\references{
  Zeileis A., Hornik K. and Murrell P. (2009),
  Escaping RGBland: Selecting Colors for Statistical Graphics.
  \emph{Computational Statistics & Data Analysis}, Forthcoming.
  Preprint available from \url{http://statmath.wu-wien.ac.at/~zeileis/papers/Zeileis+Hornik+Murrell-2008.pdf}.
}

\seealso{\code{\link{paircomp}}}

\examples{
data("CEMSChoice")
par(mar = c(5, 6, 3, 6))
plot(CEMSChoice$preference, abbreviate = FALSE)
}

\keyword{classes}
