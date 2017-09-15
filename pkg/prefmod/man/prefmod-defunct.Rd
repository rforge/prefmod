\encoding{UTF-8}
\name{prefmod-defunct}
\alias{prefmod-defunct}
\alias{plotworth}
\title{Defunct Functions in Package \pkg{prefmod}}
\description{%
  A list of functions that are no longer part of \pkg{prefmod}.
}
\usage{
plotworth(worthmat, main = "Preferences", ylab = "Estimate", psymb = NULL,
        pcol = NULL, ylim = range(worthmat), \dots)
}
\arguments{%
  \item{worthmat}{parameter matrix as generated from \code{\link{llbt.worth}} or \code{\link{patt.worth}}.}
  \item{main}{main title of the plot.}
  \item{ylab}{y-axis label}
  \item{psymb}{plotsymbols for objects, see Details below}
  \item{pcol}{colours for objects, see Details below}
  \item{ylim}{limits for y-axis}
  \item{\dots}{further graphical parameters, use e.g. \code{log = "y"} to obtain a logarithmic plot}
}
\details{%
  \code{\link{plotworth}()} was initially used to plot worth or model parameters from LLBT or pattern models (in a matrix created by \code{\link{llbt.worth}()} or \code{\link{patt.worth}()}).
  Now, the generic \code{\link[=plot.wmat]{plot}} (i.e., \code{\link{plot.wmat}()}) has to be used.
}
%\seealso{%
%}
\examples{
\dontrun{
###################
### plotworth() ###
###################

# fit only first three objects with SEX effect
m2 <- pattPC.fit(cemspc, nitems = 3, formel = ~SEX, elim = ~SEX, undec = TRUE)

# calculate and plot worth parameters
m2worth <- patt.worth(m2)
plot.wmat(m2worth)
plot.wmat(m2worth, pcol = "terrain")

# display of some plotsymbols (pch)
plot(0:25, rep(1, 26), pch = 0:25, cex = 1.5)
text(0:25, rep(0.95, 26), 0:25)

# usage of the "RColorBrewer" package
library("RColorBrewer")
mypalette <- brewer.pal(3, "Set1")
plot.wmat(m2worth, pcol = mypalette)}
}
