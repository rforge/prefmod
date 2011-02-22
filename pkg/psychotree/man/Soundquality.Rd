\name{Soundquality}
\alias{Soundquality}

\title{Quality of Multichannel Reproduced Sound}

\description{
  Paired comparison judgments of 40 selected listeners with respect to
  eight audio reproduction modes and for types of music. 
}

\usage{data("Soundquality")}

\format{
A data frame containing 783 observations on 6 variables.
\describe{
  \item{id}{Factor. Listener ID.}
  \item{time}{Factor. Listening experiment before or after elicitation and
    scaling of more specific auditory attributes.}
  \item{progmat}{Factor. The Program material: Beethoven, Rachmaninov, Steely
    Dan, Sting.}
  \item{repet}{The repetition within each time point.}
  \item{session}{The experimental session coding the presentation order.}
  \item{preference}{Paired comparison of class \code{\link{paircomp}}.
    Preferences for all 28 paired comparisons from 8 audio reproduction modes:
    Mono, Phantom Mono, Stereo, Wide-Angle Stereo, 4-channel Matrix,
    5-channel Upmix 1, 5-channel Upmix 2, and 5-channel Original.}
}
}

\details{
  See Choisel and Wickelmaier (2007) for details.
}

\references{
  Choisel, S. & Wickelmaier, F. (2007). Evaluation of multichannel
  reproduced sound: Scaling auditory attributes underlying listener
  preference. Journal of the Acoustical Society of America, 121(1),
  388-400. 
}

\seealso{\code{\link{paircomp}}}

\examples{
data("Soundquality", package = "psychotree")
summary(Soundquality$preference)
ftable(xtabs(~ time + repet + progmat, data = Soundquality))
}

\keyword{datasets}
