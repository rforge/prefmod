########################################################################
# llbt.design.R                                                        #
#                                                                      #
# version 0.1   15.2.2008  based on OBTdesign051                       #
#                          input now through function call             #
#                          options for ctrl object, dataframe, filename#
#                          default values in call (to be overwritten   #
#                             either by pasing parameters or ctrl obj) #
#                           ctrl object unified with pcpatt0           #
#                          only REDUCECAT !!                           #
#                                                                      #
# ---------------------------------------------------------------------#

llbt.design<-function(obj, nitems=NULL, objnames="", blnCasewise=FALSE, cov.sel="",
       blnGLIMcmds=FALSE,glimCmdFile="",outFile="")
{

#
# initialisations
#
  dfr<-NULL
  if(is.character(obj)){
          datafile<-obj
          nobj         <- nitems
          objnames     <- objnames
          casewise     <- blnCasewise
          glimoutput   <- blnGLIMcmds
          glimoutFile  <- glimCmdFile
          outFile      <- outFile
  } else if(is.data.frame(obj)){
          dfr<-obj
          nobj         <- nitems
          objnames     <- objnames
          casewise     <- blnCasewise
          glimoutput   <- blnGLIMcmds
          glimoutFile  <- glimCmdFile
          outFile      <- outFile
  } else if (is.list(obj)) {
          datafile     <- obj$datafile
          nobj         <- obj$nitems
          objnames     <- obj$objnames
          #reduce.cat   <- obj$blnReducecat
          casewise     <- obj$blnCasewise
          cov.sel      <- obj$cov.sel
          glimoutput   <- obj$blnGLIMcmds
          glimoutFile  <- obj$glimCmdFile
          outFile      <- obj$outFile
  } else {
          stop("first argument must be either ctrlobject, datafilename or dataframe")
  }

if (is.null(nobj)) stop("number of items not defined")



    reduce.cat <- TRUE # always TRUE until clarification

    if (length(objnames)!=nitems)
       objnames     <- paste("o",1:nobj,sep="")
    ncomp <- nobj * (nobj-1) / 2


###################only individual data for the time being
###if (aggreg) {                       # aggregated data
###
###   y<-scan(datafile)
###   if (reduce.cat) y <- red.cat.agg(y,nrespcat.data)
###
###} else {                            # individual data
####################################################

    if (is.null(dfr)) {
        if(file.access(datafile, mode=0) == 0){
           dat<-as.matrix(read.table(datafile,header=TRUE))  # datafile
        } else {
           stop("\ninput data file does not exist!\n")
        }
    } else {
        dat<-as.matrix(dfr)                                  # dataframe
        dat<-apply(dat,2,as.numeric)
    }
    nrespcat<-diff(range(dat[,1:ncomp],na.rm=TRUE))   # number of response categories categories
    raw.data <- dat[,1:ncomp]


# reduce categories
    if (min(raw.data,na.rm=TRUE)<0) raw.data<- -raw.data
    # data shifted to 0,...
    raw.data <- raw.data - min(raw.data,na.rm=TRUE)

    # check for cases with -1 1
    tab<-table(raw.data)
    levs<-as.numeric(names(tab))


      nlevs<-length(levs)
      mid<-median(levs)

       if(nlevs%%2==0) {                           # no undecided

             rawdata<-ifelse(raw.data<mid,0,1)
             nrespcat<-1
       } else {
             rawdata<-ifelse(raw.data<mid,0,ifelse(raw.data>mid,2,1))
             nrespcat<-2
       }



    nsubj<-nrow(dat)
    nrows <- ncomp * (nrespcat+1)        # number of rows for 1 design matrix (1 subject or no covs)

    blnSubjcov<-cov.sel[1]!=""             #    FALSE if ""

    if (blnSubjcov) {
       inpcovnames <- colnames(dat)[(ncomp+1):ncol(dat)]
       if (toupper(cov.sel[1]) == "ALL"){   # all covariates included
          cov.sel <- inpcovnames
       } else if(length(setdiff(cov.sel,inpcovnames))>0) {
           stop("\ncovariate name(s) in cov.sel incorrectly specified\n")
       }
       cov.case<-as.matrix(dat[,c(cov.sel)])
       #if (any(is.na(cov.case)))
       #    stop("subject covariates with NAs not allowed")

       NAs<-which(!complete.cases(cov.case))                  # check for NA
       if (length(NAs)>0){

           warning("\tselected subject covariates have NAs in lines ",paste(NAs,collapse=",")," - removed from data\n")
           notNAs<-which(complete.cases(cov.case))
           dat<-dat[notNAs,]
           cov.case<-as.matrix(cov.case[notNAs,])
           nsubj<-nrow(dat)
       }


       colnames(cov.case)<-cov.sel
       covnames<-cov.sel
       covlevels<-apply(cov.case,2,max)
       ncov<-length(cov.sel)
       if (any(apply(cov.case,2,min)<1))
           warning("subject covariates with values < 1, if these are factors recode them")
    } else {
       cov.case=NULL
       covlevels=NULL
       covnames=NULL
       ncov<-0
    }


if (reduce.cat){                            #?????
   nrespcat <- 1 + (nrespcat+1)%%2
   nrows <- ncomp * (nrespcat+1)
}

# initialisations for covariates

if (ncov == 0 ) {         # no subject covariates
      totlev <- 1
      ones.totlev <- 1
      if (casewise){     # metric (and categorical subject covs)
        totlev<-nsubj
        ones.totlev <- rep(1,nsubj)
      }
} else if (casewise){     # metric (and categorical subject covs)
      totlev<-nsubj
      ones.totlev <- rep(1,nsubj)

} else {                  # categorical subject covariates
      totlev <- prod(covlevels)  # total number of covariate levels
      ones.totlev<-rep(1,totlev) # vector for kronecker products

      indx <- ncov:1
      if (ncov == 1) {
           baslev <- nrows
      } else {
           baslev <- c(covlevels[2:ncov],nrows)
      }
      levmult <- rev(cumprod(baslev[indx]))

# generate subject covariates

      cov<-NULL
      for (j in 1:ncov) {
          scov<-gl(covlevels[j],levmult[j],totlev*nrows)
          cov<-cbind(cov,scov)
      }
}

ones.ncomp<-rep(1,ncomp)         # vector for kronecker products
ones.nrows<-rep(1,nrows)         # vector for kronecker products



#
# design matrix for objects
#




obj<- objdesign(nrows,nobj,nrespcat)# /nrespcat  ## to make design -1 0 1 instead of -2 0 2 in case of undecided
obj<- ones.totlev %x% obj        # stack object design matrix totlev times


#
# mu - factor for comparisons
#
mu<-rep.int(1:ncomp, rep.int((nrespcat+1),ncomp))
mu<-factor(rep(mu,totlev))           # stack mu totlev times


#
# design matrix for gammas (category parameters)
#

g <- ones.ncomp %x% diag(nrespcat+1)  # stack gamma matrix ncomp times for 1 subject
#print(g)
g <- ones.totlev %x% g              # stack gammas totlev times
#print(g)

gamnames<-paste("g", 0:nrespcat, sep = "")



# case: no subject covariates
#
if (ncov == 0 && !casewise) {

      y<-rep(0,nrows)
      for (i in 1:nsubj) {
        k<-1
        for (j in 1:ncomp) {
           t<-rawdata[i,j]
           if (!is.na(t)) {
              y[k+t]=y[k+t]+1
           }
           k <- k + nrespcat + 1
        }
      }


# case: categorical subject covariates
#
} else if (ncov>0 && !casewise) {

      y<-rep(0,totlev * nrows)
      for (i in 1:nsubj) {
          y.address <- sum((cov.case[i,]-1)*levmult)
          for (j in 1:ncomp) {
             t<-rawdata[i,j]
           if (!is.na(t)) {
                y[y.address+t+1]=y[y.address+t+1]+1
             }
          y.address <- y.address + nrespcat+1
          }
      }

# case: metric (and categorical) subject covariates
#
} else {

   if(ncov==0){
      cov<-NULL
   } else {
      cov<-cov.case %x% ones.nrows  # expand covariates for all response categories
      cov<-data.frame(cov)
   }
      case<-rep.int(1:nsubj, rep.int(nrows,nsubj))
      case<-factor(case)                            # factor for cases
      cov<-cbind(cov,case)
      cov<-data.frame(cov)
      cov$case<-factor(cov$case)
      covnames<-c(covnames,"CASE")
      covlevels<-c(covlevels,nsubj)
      k<-1
      y<-rep(0,nrows * totlev)

      for (i in 1:nsubj) {
        for (j in 1:ncomp) {
           t<-rawdata[i,j]
           if (!is.na(t)) {
              y[k+t]=y[k+t]+1
           }
           k <- k + nrespcat + 1
        }
      }

}


#
# prepare dataframe and export
#

if (nrespcat%%2==0) obj<-obj/2   ### just as long as ordinal model not clarified
                                 ### converts -2 0 2 to -1 0 1
if (ncov == 0 && !casewise) {
      dm<-data.frame(y,mu,g,obj)
      varnames<-c("mu",gamnames,objnames)
} else {
      dm<-data.frame(y,mu,g,obj,cov)
      varnames<-c("mu",gamnames,objnames,covnames)
}



################################################################
#
# GLIM commands output
#

if (glimoutput) {

    names(dm[,1])<-"!y"
    write.table(dm,outFile,quote=F,row.names=F)  #

    txt<-paste("$SL ",length(y),sep="")
    write(txt,file=glimoutFile)

    txt<-paste("$DATA y ",paste(varnames,collapse=" "),sep="")
    write(txt,file=glimoutFile,append=T)

    txt<-paste("$DINP '",outFile,"' 132",sep="")
    write(txt,file=glimoutFile,append=T)

    txt<-paste("$FAC mu ",ncomp," ",sep="")
    if (ncov>0) {
         facs<-na.omit(data.frame(covnames,covlevels))
         txt<-paste(txt, paste(facs[,1],facs[,2],sep = " ",collapse=" "),sep="")
    }
    write(txt,file=glimoutFile,append=T)

    if (ncov==0) {
         txt<-"$ELIM mu"
    } else if(casewise) {
         txt<-"$ELIM mu*case"
    } else {
         txt<-paste("$ELIM mu*",paste(covnames,collapse="*"),sep="")
    }
    write(txt,file=glimoutFile,append=T)

    write("$ERR P $YV y",file=glimoutFile,append=T)

    txt<-paste("$FIT",paste(objnames,collapse="+"),sep=" ")
    write(txt,file=glimoutFile,append=T)

    write("$DISP E",file=glimoutFile,append=T)

    write("$RETURN",file=glimoutFile,append=T)
}

names(dm)<-c("y",varnames)
if (glimoutput) {
  invisible(dm)
} else {
  return(dm)
}

}
