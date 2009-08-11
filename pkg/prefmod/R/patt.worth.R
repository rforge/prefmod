patt.worth<-function(obj, obj.names=NULL, outmat="worth")
{
    if(class(obj)!="pattMod") stop("function only for objects of class pattMod (see help e.g. for pattPC.fit)")

    envList<-obj$envList
    ncovpar<-envList$ncovpar

    Tmod<- regexpr("T",envList$resptype)>0               # check if time model

    if(!Tmod){                                           # if not a time model
         tpoints<-1
         nobj<-obj$envList$nobj
         npar<-(nobj - 1) * ncovpar
         lambda<-obj$coefficients[1:npar]
         lmat<-matrix(lambda,nr=nobj-1)
         lmat<-rbind(lmat,rep(0,ncol(lmat)))
    } else {                                             # time model
         tpoints<-envList$tpoints
         nobj<-(envList$nitems-1)
         npar<-nobj*ncovpar*tpoints
         lambda<-obj$coefficients[1:npar]
         lmat<-matrix(lambda,nr=nobj)
         lmat<-rbind(lmat,rep(0,ncol(lmat)))
         lmat<-matrix(lmat,nc=tpoints)
         nobj<-envList$nitems*tpoints
    }

    covlevels<-obj$envList$covlevels
    if (is.null(covlevels)) covlevels<-1

    ## preference parameters (summed lambdas) for cov groups
    struct <- unique(obj$envList$covdesmat)

    if (ncol(struct)==0) struct<-as.matrix(1)
    # summation matrix for objects
    dd<-diag(nobj)
    sum.mat <- struct %x% dd
    # sum up
    group.est <- sum.mat %*% as.vector(lmat)

    ## labels for cov groups
    if(!Tmod){                                           # if not a time model
         x<-obj$envList$model.covs
         if (is.null(x)) {
              gr.labels <- "estimate"
         } else{
         xx<-mapply(function(x,y)paste(x,y,sep=""), colnames(x),data.frame(x))
         gr.labels <-apply(xx,1,paste,collapse=":")
         }
    } else
         gr.labels <- ""




    mltp<-2
    worthmatrix<-NULL
    est<-matrix(group.est,nr=nobj/tpoints)

    ## worth matrix
    for (i in 1:ncol(est)) {

       # worth parameters

       worth<-rep(0,nobj/tpoints)
       coeff<-est[,i]
       worthdenominator<-0
       for (j in 1:(nobj/tpoints)) {
         worthdenominator<-worthdenominator+exp(mltp*coeff[j])
       }
       for (j in 1:(nobj/tpoints)) {
         worth[j]<-exp(mltp*coeff[j])/worthdenominator
       }
       worthmatrix<-cbind(worthmatrix,worth)
    }

    if (is.null(obj.names)){
       obj.names<-obj$envList$obj.names[1:(nobj/tpoints)] # default: only names for first time point are used
    } else {
       obj.names<-obj.names[1:(nobj/tpoints)]
    }
    ## label worth matrix
    if(Tmod) {
       if(gr.labels[1]=="") Tlabel<-"T" else Tlabel=":T"
       worth.names<-paste(rep(gr.labels,rep(tpoints,length(gr.labels))),paste(Tlabel,1:tpoints,sep=""),sep="")
       colnames(worthmatrix)<-worth.names
       rownames(worthmatrix)<-obj.names
    } else {
       colnames(worthmatrix)<-gr.labels
       rownames(worthmatrix)<-obj.names
    }

    colnames(est) <- colnames(worthmatrix)
    rownames(est) <- rownames(worthmatrix)

    #class(worthmatrix) <- c("pattW")                         #class: pattern worth
    #worthmatrix

    switch(outmat,
       "lambda" = return(est),
       "worth" = return(worthmatrix),
       # "est" = return(lambda.mat),
       stop("     outmat must be either 'worth' or 'lambda'\n")
    )

}
