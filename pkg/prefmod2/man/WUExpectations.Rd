\name{WUExpectations}
\alias{WUExpectations}
\encoding{latin1}

\title{Expectations about Assistants and Professors at WU Wien}

\description{
  Data from a survey of 96 assistants and 31 professors at WU Wien
  examining which skills of assistants and professors are preferred.
}

\usage{data("WUExpectations")}

\format{
A data frame containing 127 observations on 7 variables.
  \describe{
    \item{\code{position}}{Factor coding position at university (either assistant or professor).}
    \item{\code{gender}}{Factor coding gender.}
    \item{\code{parttime}}{Factor. Did the interviewee work parttime or fulltime?}
    \item{\code{seniority}}{Numeric. Work experience at university in years.}
    \item{\code{deptsize}}{Numeric. Number of employees within the department.}
    \item{\code{professor}}{Paired comparison of class \code{\link[psychotree]{paircomp}}.
      Preferences for all 28 paired comparisons from 8 professor skills: competence,
      trustworthiness, communication, fairness, reliability, interest, motivation,
      cooperation.}
    \item{\code{assistant}}{Paired comparison of class \code{\link[psychotree]{paircomp}}.
      Preferences for all 28 paired comparisons from 8 assistant skills: competence,
      reliability, interest, creativity, cooperation, trustworthiness, diligence,
      motivation.}
  }
}

\details{
Pongratz and Weiler (2007) carried out the survey about
preferences on skills of assistants and professors at WU Wirtschaftsuniversität Wien
(Vienna University of Economics and Business). During a presurvey at Universität Innsbruck
the eight main important personal skills for being a professor and for being an
assistant, respectively, had been determined.

Professors and assistants at the WU were presented with all 28 pairs of
professor skills and with all 28 pairs of assistant skills. Each pair (AB) had
to be judged on the following scale: A is much more important than B; A is more
important than B; A and B are equal; B is more important than A; B is much
more important than A. Additionally to the preferences concering the skills
questions regarding gender, number of employees within the department and work
experience at a university in years were asked.

The data were reanalyzed by Weber and Weiler (2009) focusing on the question
whether the expectations about professors (or assistants) coinced or differed
with respect to the interviewees own position. A similar analysis is included in
the example section below. Unfortunately, the results can only be reproduced approximately
as their preprocessing is not fully reconstructable. However, \code{WUExpectations} contains
the full and correct version of the data. The preprocessed data as used by Weber and
Weiler (2009) is available on the accompanying URL together with their replication
code.
}

\references{
Pongratz C., Weiler M. (2007).
  \emph{Der Paarvergleich als Skalierungsmethode zur Messung von wechselseitigen
  Erwartungen zwischen Universitätsbediensteten.}
  Diploma Thesis, WU Wirtschaftsuniversität Wien.

Weber D., Weiler M. (2009).
  Erwartungskoinzidenz bei UniversitätsprofessorInnen und wissenschaftlichen
  MitarbeiterInnen: Das log-lineare Bradley-Terry Modell für Paarvergleichsdaten.
  In R. Hatzinger, R. Dittrich, T. Salzberger (eds.),
  \emph{Präferenzanalyse mit R}, Facultas Universitätsverlag, Wien.
  URL \url{http://prefmod.R-Forge.R-project.org/PAmitR/}
}

\author{Daniela Weber}

\seealso{
 \code{\link{btreg}, \code{\link[psychotree]{paircomp}}}
}

\examples{
## data and preprocessing
data("WUExpectations", package = "prefmod2")

## subset selection: use only complete cases
wu <- na.omit(WUExpectations)

## recode comparisons from likert scale
mscale(wu$assistant) <- c(-1, -1, 0, 1, 1)
mscale(wu$professor) <- c(-1, -1, 0, 1, 1)



## simple Bradley-Terry regression models

## expectations about assistants
asnt_by_asnt <- btreg(assistant ~ 1, data = wu, subset = position == "assistant")
asnt_by_prof <- btreg(assistant ~ 1, data = wu, subset = position == "professor")
summary(asnt_by_asnt)
summary(asnt_by_prof)

## expectations about professors
prof_by_asnt <- btreg(professor ~ 1, data = wu, subset = position == "assistant")
prof_by_prof <- btreg(professor ~ 1, data = wu, subset = position == "professor")
summary(prof_by_asnt)
summary(prof_by_prof)



## visualization via parallel coordinates
par(mfrow = c(1, 2), mar = c(3.5, 2, 3, 0.5))
plot(asnt_by_asnt, main = "Expectations about assistants",
  abbreviate = 4, ylim = c(0, 0.3), ylab = "")
lines(worth(asnt_by_prof), lty = 3, type = "b", pch = 2, col = 2)
legend("topright", c("by assistants", "by professors"),
  lty = 3:2, pch = 1:2, col = 1:2, bty = "n")
plot(prof_by_asnt, main = "Expectations about professors",
  abbreviate = 4, ylim = c(0, 0.3), ylab = "")
lines(worth(prof_by_prof), lty = 3, type = "b", pch = 2, col = 2)
legend("topright", c("by assistants", "by professors"),
  lty = 3:2, pch = 1:2, col = 1:2, bty = "n")


## visualization via scatterplots
par(mfrow = c(1, 2), mar = c(4.5, 4, 4, 0.5))
plot(worth(asnt_by_asnt), worth(asnt_by_prof),
  xlim = c(0, 0.3), ylim = c(0, 0.3), type = "n",
  main = "Expectations about assistants",
  xlab = "by assistants", ylab = "by professors")
text(worth(asnt_by_asnt), worth(asnt_by_prof), names(worth(asnt_by_asnt)))
abline(0, 1)
plot(worth(prof_by_asnt), worth(prof_by_prof),
  xlim = c(0, 0.3), ylim = c(0, 0.3), type = "n",
  main = "Expectations about professors",
  xlab = "by assistants", ylab = "by professors")
abline(0, 1)
text(worth(prof_by_asnt), worth(prof_by_prof), names(worth(prof_by_asnt)))
}

\keyword{datasets}
