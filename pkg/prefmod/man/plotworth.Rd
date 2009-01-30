\name{plotworth}
\alias{plotworth}
\title{Function to plot worth parameters from pattern models}
\description{
  A plot of the worth parameter matrix obtained from the fit of a pattern model is produced.
  The worth parameter matrix is produced by \code{\link{patt.worth}}.
}
\usage{
plotworth(worthmat, main = "Preferences", ylab = "Estimate", psymb = NULL, ylim = range(worthmat))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{worthmat}{worth parameter matrix, object of class \code{pattW}.
       Alternatively, the contents of a user defined matrix can also be plotted
       (for an example see \code{\link{llbt.design}},
       provided that the rows correspond to objects and the columns to groups.}
  \item{main}{main title of the plot.}
  \item{ylab}{y-axis label}
  \item{psymb}{plotsymbols for objects, see Details below}
  \item{ylim}{limits for y-axis}
}
\details{
  Plotsymbols can be defined as an integer vector of length equal to the
  number of objects, e.g., \code{psym = c(15,22,18)}. They specify the graphical option \code{pch} as used in
  the \code{\link{points}} function. The default (\code{psymb = NULL}) uses
  the symbols 15 through 18 and 21 through 25. The number of symbols is determined
  from the number of rows in \code{worthmat}. A display of some plotsymbols may
  be obtained from the correponding example below.
 }
\author{Reinhold Hatzinger}
\seealso{\code{\link{patt.worth}}}
\examples{
    ## display of some plotsymbols (pch)
    plot(0:25,rep(1,26),pch=0:25,cex=1.5)
    text(0:25,rep(0.95,26),0:25)

    ## fit only first three objects with SEX effect
    data(cemspc)
    m2<-pattPC.fit(cemspc, nitems=3, formel=~SEX, elim=~SEX, undec=TRUE)

    ## calculate and plot worth parameters
    m2worth<-patt.worth(m2)
    plotworth(m2worth)
    }
\keyword{models}
