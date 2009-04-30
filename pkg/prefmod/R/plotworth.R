plotworth<-function(worthmat, main="Preferences", ylab="Estimate",
              psymb=NULL, pcol=NULL, ylim = range(worthmat))
{
#
# plot ranking
#
#


coeff<-unclass(worthmat)
objnames<-rownames(coeff)
grnames<-colnames(coeff)

# for proper spacing x-axis labels (labels are broken into lines at ":")

if(length(grep(":",grnames))){
  gr.lines<-unlist(strsplit(grnames[1],":"))
  ngrlin<-length(gr.lines)-0.5
  grnames<-gsub(":","\n",grnames)
} else ngrlin <- 1

nobj<-dim(coeff)[1]
ngroups<-dim(coeff)[2]
if (ngroups == 1) colnames(coeff) <- ""

# plotsymbols and color

if (is.null(psymb)) psymb<-c(15:18,21:25)[1:nobj]

if (is.null(pcol))
     farbe <- rainbow(nobj)
else if (length(pcol)>1)
     farbe <- pcol
else if(pcol=="black")
     farbe <- "black"
else if(pcol %in% c("heat","topo","terrain","cm","gray"))
     farbe <- eval(call(paste(pcol,".colors",sep="",collapse=""), nobj))
else farbe <- rainbow(nobj)

#     gray(0:(nobj-1) /  nobj)


#  farbe<-rainbow(26, s = 1, v = 1, start = 0, end = 1, gamma = 1)
#  farbe <- brewer.pal(nclass = 26,  palette = "Spectral")
#  plot(0:25,rep(1,26),pch=0:25,col=farbe,cex=1.5)
#  text(0:25,rep(0.95,26),0:25)
#

## plot
par(omi = c(0.2,0.2,0.5,0.2), mar=c(3,4,0.1,0) ) # makes plot fill the whole region

plot(c(0.5,ngroups+0.5),c(min(coeff),max(coeff)),type="n",axes=FALSE, xlab="",ylab=ylab,ylim=ylim)
title(main,outer=TRUE)
box()
axis(2)
axis(1,at=1:ngroups,labels=grnames,mgp=c(3,ngrlin,0))

adj<-1/nchar(as.character(objnames))
adj<-strwidth(as.character(objnames))/2
d<-strwidth("d")                        # distance from point to label
adj<-adj+d


for (i in 1:ngroups) {
  pm<-rep(c(1,-1),nobj)
  pm<-pm[1:nobj]                 # plus/minus orientation of point labels
  o<-order(coeff[1:nobj,i])
  sadj<-adj[o]
  sobjnames<-objnames[o]         # sorted obj names
  spsymb<-psymb[o]               # sorted plotsymbols
  scoeff<-sort(coeff[,i])        # sorted estimates

  x<-rep(i,nobj)
  #plot(x,scoeff,type="o")
  xy <- xy.coords(x, scoeff)

  lines(c(i,i), range(scoeff), col = "gray", lty = "dotted")
  points(xy, pch=spsymb, cex=1.5, col=farbe)
  text(x+sadj*pm,scoeff,sobjnames)

}
}
