\encoding{UTF-8}
\name{music}
\alias{music}
\docType{data}
\title{Data (ratings): Music (US General social survey 1993)}
\description{%
  The dataset contains data on 18 items referring to the liking of musical styles (ratings on a 5-point Likert type response scale) and three subject covariates.
  The data is an excerpt from the US General Social Survey 1993.
  Missing values have been removed from the subject variables.
}
\usage{music}
\format{
  A data frame with 1597 observations.
  The variables \code{bigb} through \code{hvym} are are responses to items asking for liking/disliking of type of music:
  Do you like it very much (1), like it (2), have mixed feelings (3), dislike it (4), dislike it very much (5).
  \describe{
    \item{\code{bigb}}{like or dislike bigband music}
    \item{\code{blug}}{like or dislike bluegrass music}
    \item{\code{coun}}{like or dislike country western music}
    \item{\code{blue}}{like or dislike blues or r and b music}
    \item{\code{musi}}{like or dislike broadway musicals}
    \item{\code{clas}}{like or dislike classical music}
    \item{\code{folk}}{like or dislike folk music}
    \item{\code{gosp}}{like or dislike gospel music}
    \item{\code{jazz}}{like or dislike jazz}
    \item{\code{lati}}{like or dislike latin music}
    \item{\code{mood}}{like or dislike easy listening music}
    \item{\code{newa}}{like or dislike new age music}
    \item{\code{oper}}{like or dislike opera}
    \item{\code{rap}}{like or dislike rap music}
    \item{\code{regg}}{like or dislike reggae music}
    \item{\code{conr}}{like or dislike contemporary rock music}
    \item{\code{oldi}}{like or dislike oldies rock music}
    \item{\code{hvym}}{like or dislike heavy metal music}
    \item{\code{age}}{age in years}
    \item{\code{educ}}{highest year of school completed}
    \item{\code{sex}}{1 male, 2 female}
  }
}
\source{%
  Dataset: General Social Survey 1993 [United States]\cr%
  CPANDA Identification Number: a00006\cr%
  Cultural Policy and the Arts National Data Archive\cr%
  \url{http://www.cpanda.org/cpanda/getDDI.xq?studyID=a00006#dataDscr}\cr%
}
\examples{
summary(music)
}
\keyword{datasets}
