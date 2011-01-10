checkMIS<-function(obj,nitems,MISmodel="obj",verbose=FALSE){


  nobj<-nitems
  ## prepare data
  ncomp<-nobj*(nobj-1)/2
  if(is.character(obj)){                    ## datafilename supplied
        datafile    <-  obj
        if(file.access(datafile, mode=0) == 0){
           dat<-as.matrix(read.table(datafile,header=TRUE))  # datafile
        } else {
           stop("\ninput data file does not exist!\n")
        }
   } else if(is.data.frame(obj)){            ## data frame supplied
        dat<-as.matrix(obj)                                  # dataframe
        dat<-apply(dat,2,as.numeric)
   } else {
        stop("first argument must be either datafilename or dataframe")
   }
   varnames<-colnames(dat)
   if (ncol(dat)>ncomp) {
        covnames<-varnames[(ncomp+1):ncol(dat)]
        covs<-as.data.frame(dat[,(ncomp+1):ncol(dat)])
   } else {
        covs<-NULL
   }

   ## vector of counts of NAs, columnwise
   misv<-apply(dat,2,function(x) sum(is.na(x)) )
   mism<-matrix(0,nobj,nobj)
   mism[upper.tri(mism,diag=FALSE)]<-misv
   mism<-mism+t(mism)
   rownames(mism)<-colnames(mism)<-paste("o",1:nobj,sep="")
   objMis<-colSums(mism)
   names(objMis)<-paste("o",1:nobj,sep="")
   if(verbose){
     cat("number of missing comparisons:\n")
     print(mism)
     cat("number of missing comparisons for objects:\n",objMis,"\n")
   }
   if (MISmodel=="obj")
      RET<-ifelse(objMis==0,F,T)
   else if (MISmodel=="comp")
      RET<-ifelse(misv==0,F,T)
   else
      stop('\nNo output: MISmodel not correctly specified. Use "obj" or "comp"\n')

   invisible(RET)
}
