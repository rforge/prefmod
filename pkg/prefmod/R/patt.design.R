`patt.design` <- local({

ENV<-new.env()  # environment local to pcpatt0 - no objects in .GlobalEnv

###function(obj,dfr=NULL)
function(obj, nitems=NULL, objnames="", resptype="paircomp", blnRevert=FALSE, cov.sel="",
       blnIntcovs=FALSE, blnGLIMcmds=FALSE,glimCmdFile="",outFile="", intFile="")
{
  env2<-new.env()

  #' datafile     = "",       # dataframe used
  #' nitems       = 4,
  #' blnRevert    = FALSE,
  #' blnReducecat = TRUE,
  #' blnIntcovs   = TRUE,
  #' resptype    = "ranking",
  #'
  #' cov.sel      = "",
  #'
  #' blnGLIMcmds  = FALSE,    # no GLIM output
  #' glimCmdFile  = "",
  #' outFile      = "",
  #' intFile      = ""

  # default ctrl object
  ctrl<-list(datafile="", nitems=NULL, objnames="", resptype="paircomp", blnRevert=FALSE, cov.sel="",
       blnIntcovs=FALSE, blnGLIMcmds=FALSE,glimCmdFile="",outFile="", intFile="")

  dfr<-NULL
  if(is.character(obj)){                    ## datafilename supplied
          ctrl$datafile    <-  obj
          ctrl$nitems      <-  nitems
          ctrl$objnames    <-  objnames
          ctrl$resptype    <-  resptype
          ctrl$blnReducecat <- TRUE         # always TRUE until clarification
          ctrl$blnRevert   <-  blnRevert
          ctrl$cov.sel     <-  cov.sel
          ctrl$blnIntcovs  <-  blnIntcovs
          ctrl$blnGLIMcmds <-  blnGLIMcmds
          ctrl$glimCmdFile <-  glimCmdFile
          ctrl$outFile     <-  outFile
          ctrl$intFile     <-  intFile
  } else if(is.data.frame(obj)){            ## data frame supplied
          dfr<-obj
          ctrl$datafile    <-  ""
          ctrl$nitems      <-  nitems
          ctrl$objnames    <-  objnames
          ctrl$resptype    <-  resptype
          ctrl$blnReducecat <- TRUE         # always TRUE until clarification
          ctrl$blnRevert   <-  blnRevert
          ctrl$cov.sel     <-  cov.sel
          ctrl$blnIntcovs  <-  blnIntcovs
          ctrl$blnGLIMcmds <-  blnGLIMcmds
          ctrl$glimCmdFile <-  glimCmdFile
          ctrl$outFile     <-  outFile
          ctrl$intFile     <-  intFile
  } else if (is.list(obj)) {                ## ctrl object
          for (i in names(obj))             #       replaces default values
                  ctrl[[i]]<-obj[[i]]
          ctrl$blnReducecat <- TRUE         # always TRUE until clarification
  } else {
          stop("first argument must be either ctrlobject, datafilename or dataframe")
  }

     ####from pcpatt0 - obsolete
     ####if (!is.list(ctrl))
     ####   stop(paste(deparse(substitute(ctrl)),"is not a list object - see help for pcpatt0"))

     ## initialising
     ##
     cat("initialising...\n")
     flush.console()

     for (i in 1:length(ctrl))
        do.call("assign",list(names(ctrl)[i],ctrl[[i]],ENV))

     init(dfr)    # initialisation function
     dat<-ENV$dat

     ## all possible response patterns and/or difference patterns
     ##
     cat("generating response patterns...\n")
     flush.console()

     resptype <- get("resptype",ENV)

     if (resptype == "rating") {
        generateLpatterns(env2)
     } else if(resptype == "ranking") {
        generateRpatterns(env2)
     } else if(resptype == "paircomp") {
        generatePCpatterns(env2)
     }
     datStr   <- env2$datStr
     dpattStr <- env2$dpattStr
     npatt    <- env2$npatt
     diffs    <- env2$diffs
     blnUndec <- env2$blnUndec

     ## designmatrix-kernel for objects, undecided/categories, interactions
     ##
     cat("setting up the design matrix...\n")
     blnIntcovs <-get("blnIntcovs",ENV)
     onedesign<-designkernel(diffs,blnUndec,blnIntcovs,env2)

     # tidy up
     rm(diffs)
     rm(dat,envir=ENV)
     gc(FALSE)

     ## complete design for categorical subject covariates
     ##
     blnSubjcov <-get("blnSubjcov",ENV)
     if (!blnSubjcov) {
         ## no subject covariates
         ##
         totlev <- 1
         ones.totlev <- 1
     } else {
         ## subject covariates
         ##
         covlevels <-get("covlevels",ENV)
         covnames  <-get("covnames",ENV)
         cov.case  <-get("cov.case",ENV)
         ncov      <-get("ncov",ENV)
         totlev <- prod(covlevels)
         # vector for kronecker products (to stack design matrix)
         ones.totlev<-rep(1,totlev)
         indx <- ncov:1
         if (ncov == 1) {                          # only 1 covariate
              baslev <- npatt
         } else {                                  # >1 covariates
              baslev <- c(covlevels[2:ncov],npatt)
         }
         levmult <- rev(cumprod(baslev[indx]))
         # transform subject covariates data into covariate vectors
         scov<-matrix(0,nrow=totlev*npatt,ncol=ncov)
         colnames(scov)<-covnames
         scov<-data.frame(scov)
         for (j in 1:ncov) {
            scov[,j]<-gl(covlevels[j],levmult[j],totlev*npatt)
         }
     }


     ## extension of design-kernel in case of subject covariates
     ##
     design<-ones.totlev %x% onedesign
     colnames(design)<-colnames(onedesign)
     # tidy up
     rm(onedesign)
     gc(FALSE)


     ## calculate response pattern frequencies based on
     ## comparison of string representation of
     ## possible patterns and observed patterns
     ##
     cat("calculating response pattern frequencies...\n")
     flush.console()

     # count occurrency of patterns into y
     # (according to covariates if blnSubjcov==TRUE)
     nsubj <-get("nsubj",ENV)
     y<-rep(0,totlev * npatt)
     cov.addr <- 0       # for case blnSubjcov==FALSE
     for (i in 1:nsubj) {
         if (blnSubjcov)
            cov.addr <- sum((cov.case[i,]-1)*levmult)
         j <- match(datStr[i],dpattStr)
         y[j+cov.addr]<-y[j+cov.addr]+1
     }
     ##alternatively:
     #j<-sapply(datStr,function(x)match(x,dpattStr),USE.NAMES=FALSE)
     #if(blnSubjcov){
     #    cov.addr<-apply(cov.case,1,function(x)crossprod((x-1),levmult))
     #    tb<-table(j+cov.addr)
     #} else {
     #    tb<-table(j)
     #}
     #y[as.numeric(names(tb))]<-tb
     #rm(j,cov.addr,tb)

     # tidy up
     rm(datStr,dpattStr)
     rm(cov.case, envir=ENV)
     gc(FALSE)


     ## prepare dataframe
     ##
     dm<-data.frame(cbind(y,design))
     if (blnSubjcov){
           dm<-data.frame(dm,scov)
           rm(scov,cov.case)
     }
     # tidy up
     rm(y,design)
     gc(FALSE)

     ## generate files for GLIM
     ##
     if (get("blnGLIMcmds",ENV)){
          nintcovs.out<-40 # max number of values/line in interaction output file
          writeGLIMcmds(dm,blnUndec,
              blnIntcovs,outFile,ncov,env2$nintpars,intFile,nintcovs.out,glimCmdFile,
              covnames,covlevels,objnames,ENV$undecnames)
         #writeGLIMcmds(dm,blnUndec,ENV)
     } else {
         return(dm) # R output only if blnGLIMcmds==FALSE
     }
     cat("Done\n\n")
}
}) # end local
