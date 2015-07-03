\encoding{UTF-8}
\name{plot.wmat}
\alias{plot.wmat}
\title{Method to plot worth or model parameters from LLBT or pattern models}
\description{
  A plot of the worth or model parameter matrix obtained from the fit of an LLBT or  pattern model is produced.
This matrix is obtained from \code{llbt.worth} or \code{patt.worth} and is an object of class \code{wmat}.
}
\usage{
\method{plot}{wmat}(x, main = "Preferences", ylab = "Estimate",
       psymb = NULL, pcol = NULL, ylim = range(worthmat),log="", ...)
}
\arguments{
  \item{x}{worth or parameter matrix as generated from \code{\link{llbt.worth}} or \code{\link{patt.worth}}.}
  \item{main}{main title of the plot.}
  \item{ylab}{y-axis label}
  \item{psymb}{plotsymbols for objects, see Details below}
  \item{pcol}{colours for objects, see Details below}
  \item{ylim}{limits for y-axis}
  \item{log}{if specified as  \code{log = "y"}, the y-axis is to be logarithmic}
  \item{...}{further graphical parameters, use e.g. \code{lty = "dashed"} to obtain dashed lines connecting the objects}
}
\details{
  Plotsymbols can be defined as an integer vector of length equal to the
  number of objects, e.g., \code{psymb = c(15,22,18)}. They specify the graphical option \code{pch} as used in
  the \code{\link{points}} function. The default (\code{psymb = NULL}) uses
  the symbols 15 through 18 and 21 through 25. The number of symbols is determined
  from the number of rows in \code{worthmat}. A display of some plotsymbols may
  be obtained from the correponding example below.

  If \code{pcol = NULL}, the colours for objects are defined from the \code{rainbow_hcl} palette using
  the \pkg{colorspace} package.
  Other specifications include "heat", "terrain" (see \code{\link[colorspace]{rainbow_hcl}}), and "gray" (see \code{\link{grDevices}}).
  The number of different colours is automatically determined via the number of objects.
  Alternatively, \code{pcol} can be specified as a character vector containing user defined RGB colour
  values for all objects (as hexadecimal strings in the form \code{"#rrggbb"}), e.g., for blue \code{"#0000FF"}).
  These are usually set up using standard colour palettes (see \code{rainbow} or, e.g., the
  \pkg{RColorBrewer} package (see Examples below).
	
	If you are still using \code{\link{plotworth}()}, please update your code!
	As of now, the routine tries to catch errors by setting \code{class(worthmat) <- "wmat"} and passing it to \code{plot.wmat(x = worthmat, ...)}.
  The function \code{plotworth()} is already considered defunct and will be removed in subsequent versions.
}
\author{Reinhold Hatzinger}
\seealso{\code{\link{patt.worth}}}
\examples{
# fit only first three objects with SEX effect
m2 <- pattPC.fit(cemspc, nitems = 3, formel = ~SEX, elim = ~SEX, undec = TRUE)

# calculate and plot worth parameters
m2worth <- patt.worth(m2)
plot(m2worth)
plot(m2worth, pcol = "terrain")

# display of some plotsymbols (pch)
plot(0:25, rep(1, 26), pch = 0:25, cex = 1.5)
text(0:25, rep(0.95, 26), 0:25)

# usage of the "RColorBrewer" package
\dontrun{
library("RColorBrewer")
mypalette <- brewer.pal(3, "Set1")
plot(m2worth, pcol = mypalette)}
}
\keyword{models}
