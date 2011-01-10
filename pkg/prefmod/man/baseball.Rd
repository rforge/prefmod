\name{baseball}
\alias{baseball}
\docType{data}
\title{Data (paired comparisons): Baseball Games}
\description{
  The result of the 1987 season for seven baseball teams in the Eastern Division of the American League according to the (home team, away team) classification are shown.
}
\usage{data(baseball)}
\format{
  Baseball is a numeric vector with the results for the season according to the (home team, away team) classification.
 }
\details{
The results of the seven teams Milwaukee, Detroit, Toronto, New York, Boston, Cleveland and Baltimore, that play 13 games each. There is no possibility of ending in a draw.
}
\references{
   Alan Agresti, \emph{Categorical Data Analysis (Second Edition)}, 2002 pages 437 and 438

   R. Dittrich, R. Hatzinger, and W. Katzenbeisser, Fitting paired comparison models in GLIM.
   \emph{GLIM newsletter} 1997

   }
\examples{
#### baseball example (Agresti, 2002, p.437)
data(baseball)

# pseudodata for generating a design matrix
d1<-c(rep(0,21),1)
d2<-c(1,rep(0,20),2)
d<-data.frame(rbind(d1,d2))
names(d)<-c(paste("v",1:21,sep=""),"cov")

# design matrix
des5<-llbt.design(d, nitems=7,
         objnames=c("MIL","DET","TOR","NY","BOS","CLE","BAL"),
         cat.scov="cov")
des5$y <- baseball
des5$mu<-gl(42,2)
pos<-c(rep(1:0, 21), rep(0:1,21) )

# fit model and display results
res5<-gnm(y~MIL+DET+TOR+NY+BOS+CLE+BAL + pos,eliminate = mu, data = des5, family = poisson)
w5<- llbt.worth(res5)
plotworth(w5)
}
\keyword{datasets}
