patt.worth<-function(obj, obj.names=NULL, outmat="worth")
{
    fitobj<-obj


#### updated on 2011-03-26
#### can now handle fitobjects from pattPCfit and
#### from gnm, provided patt.design (also new version) was used
#### object specific covariates are possible

# check if either model results from gnm (using patt.design)
# or from pattPC.fit

pattdesignmodel <- FALSE
pattPCfitmodel <- ("pattMod" %in% class(fitobj))
if (!pattPCfitmodel){
   dterm <- deparse(fitobj$call$data)            ## unten nocheinmal
   x <- tryCatch(get(dterm),error = function(e) FALSE)
   if("pattdes" %in% class(x)) pattdesignmodel <- TRUE
}
if (pattdesignmodel == pattPCfitmodel)
     stop(paste("Model result must be either from patt*.fit or having used patt.design"))


## fitobj from pattdesign, gnm
if (pattdesignmodel){

      ## subject covariates

      num.scovs<-attr(get(dterm),"num.scovs")
      if(!is.null(num.scovs))
          warning("Numerical subject covariates not (yet) implemented, they are ignored.\n
          Result is a matrix of intercepts for regression on numerical subject covariates!!")

      # extract subject covariates from design frame attribute
      cat.scovs<-attr(get(dterm),"cat.scovs")
      if(!is.null(cat.scovs))
        eterms <- cat.scovs
      else
        eterms <- vector(,0)
      # extract subject covariates from eliminate
      #eterm<-deparse(fitobj$call$eliminate)
      #if (eterm=="NULL") eterm<-"mu"         # if no eliminate specified insert a dummy mu
      #eterms<-unlist(strsplit(eterm,":"))
      #if("mu" %in% eterms)
      #   eterms<-eterms[-which(eterms=="mu")]

      ### achtung laenge 0
      if(length(eterms>0)) {

          # extract subject covariates from formula
          fterm<-deparse(fitobj$call$formula)
          fterms<-unique(unlist(strsplit(fterm,"[ ()+*:~]")))
          subjcov.names<-intersect(eterms,fterms)
          ### achtung laenge 0
          if(length(subjcov.names)>0) { ### achtung laenge 0
            # set up subj cov design
            levlist<-lapply(fitobj$model[subjcov.names],levels)
            maxlev<-sapply(levlist,length)
            subjdes<-gfac2(maxlev)
            colnames(subjdes)<-names(maxlev)
          } else {
            maxlev <- 1
            length(eterms)<-0
          }
      } else
          maxlev<-1

      ## objects

      # extract design data frame
      dterm<-deparse(fitobj$call$data)

      objs<-attr(get(dterm),"objnames")  # object names


      ## set up summation matrix

      nrep <- prod(maxlev)

      # objects
      nobj<- length(objs)
      omat <- diag(nobj)

      o<-rep(1,nrep)%x%omat
      colnames(o) <-  objs

      # summation frame
      S <- data.frame(o)

      # subjects
      if(length(eterms>0)){
          s<-subjdes %x% rep(1,nobj)
          s<-data.frame(s)
          s<-as.data.frame(lapply(s,as.factor))
          names(s)<-names(maxlev)

          # summation frame
          S <- data.frame(S,s)
      }

      # obj covs
      objcovs<-attr(get(dterm),"objcovs")
      if(!is.null(objcovs)){
         # add items to objcovs (necessary to have single objects in objcovs if required from model formula)
         ocovs<-diag(nobj)
         colnames(ocovs)<-objs
         objcovs<-cbind(ocovs,objcovs)
         # extract those obj covs which are in model formula
         vars<-attr(attr(fitobj$terms,"factors"),"dimnames")[[1]]
         vars2<-colnames(objcovs)
         diffs<-intersect(vars,vars2)
         objcovs<-objcovs[,diffs,drop=FALSE]
         oc<-rep(1,nrep)%x%objcovs
         colnames(oc)<-diffs

         # summation frame
         S <- data.frame(S,oc)
      }

      # generate summation matrix from pseudo gnm
      elim<-as.factor(1:(nobj*nrep))
      frml<-eval(fitobj$call$formula) # if formula is a symbol first eval it
      frml<-update.formula(frml, .~.-pos)     # remove position variable 2011-01-03
      ia<-names(coef(fitobj))[grep("I[0-9]+[.][0-9]+",names(coef(fitobj)))] ## remove ia variables 2011-03-30
      if (length(ia)>0) {
        remia<-paste(".~.-",paste(ia,sep="",collapse="-"),collapse="")
        frml<-update.formula(frml, remia)
      }
      u<-names(coef(fitobj))[grep("^u[0-9]+$",names(coef(fitobj)))] ## remove undec variables 2011-03-30
      if (length(u)>0) {
        remu<-paste(".~.-",paste(u,sep="",collapse="-"),collapse="")
        frml<-update.formula(frml, remu)
      }

      y<-rnorm(nrow(S)) # random y - y required for pseudo fit
      S<-data.frame(y,S,elim)

      # check which vars are unnecessary in frml
      vars<-attr(attr(fitobj$terms,"factors"),"dimnames")[[1]]
      if (any(vars=="mu")) vars<-vars[-which(vars=="mu")]
      ###if (any(vars=="pos")) vars<-vars[-which(vars=="pos")]        ####### 2011-01-03
      vars2<-colnames(S)
      diffs<-setdiff(vars,vars2)
      # add unnecessary variables to S
      if (length(diffs)>0){
         oldnc<-ncol(S)
         for (i in 1:length(diffs)) S<-data.frame(S,rnorm(nrow(S)))
         colnames(S)[(oldnc+1):(oldnc+length(diffs))]<-diffs
      }

      ## alternative: update.formula

      oldopt<-options("warn"=-1)
      sumMat<-as.matrix(model.matrix(gnm(formula=frml,eliminate=elim,data=S,iterStart = 1, iterMax = 1)))
      options(oldopt)

      coefs<-coef(fitobj)
      if (any(names(coefs)=="pos"))
            coefs<-coefs[-which(names(coefs)=="pos")] ## remove pos variable 2011-01-03
      iapos <- grep("I[0-9]+[.][0-9]+",names(coefs))  ## remove ia variables 2011-03-30
      if (length(iapos)>0) coefs <- coefs[-iapos]
      upos <- grep("^u[0-9]+$",names(coefs))          ## remove undecided variables 2011-03-30
      if (length(upos)>0) coefs <- coefs[-upos]
      if (any(names(coefs)=="(Intercept)"))
            coefs<-coefs[-1]                          ## remove 2011-03-30

      # remove NA estimates from coefficients
      notna<-!is.na(coefs)
      if (length(notna)>0){
          sumMat<-sumMat[,notna,drop=F]
          coefs<-coefs[notna]
      }

      # remove unnecessary terms from estimates and from coefficients
      if (length(diffs)>0){
        dc<-which(names(coefs) %in% diffs)
        if (length(dc)>0){    # ignore variables which are not coefficients, eg., if y is renamed 2011-01-03
          sumMat<-sumMat[,-dc, drop=FALSE]
          coefs<-coefs[-dc]
        }
      }

      sumvec<- sumMat %*% coefs

      ## lambda matrix
      lambda.groups.mat<-matrix(sumvec, nrow=nobj)


      if(length(eterms>0)) {
                 x<-subjdes


                 ## labels for cov groups
                 xx<-mapply(function(x,y)paste(x,y,sep=""), colnames(x),data.frame(x))
                 gr.labels <-apply(xx,1,paste,collapse=":")
                 colnames(lambda.groups.mat) <- gr.labels
      } else
                 colnames(lambda.groups.mat) <- "estimate"

      ## obsolete since option obj.names removed 2011-01-03
      #if(is.null(obj.names))
      #   rownames(lambda.groups.mat)<-objs
      #else
      #   rownames(lambda.groups.mat)<-obj.names
      rownames(lambda.groups.mat)<-objs


      ## in case of objcovs the rows are collapsed according to obj covs
      if (!is.null(objcovs)) {
          # collapse rownames of lambda.groups.mat
          u<-1:nobj
          names(u)<- as.numeric(factor(apply(lambda.groups.mat,1,paste,collapse="")))
          nam<-aggregate(objs,list(u[names(u)]),paste)$x

          leg<-aggregate(rownames(lambda.groups.mat),as.data.frame(objcovs),paste)
          lambda.groups.mat<-unique(lambda.groups.mat)

          if(is.matrix(nam)){
            rnam<-apply(nam,1,paste,sep="", collapse=",")
            if (nrow(lambda.groups.mat)==length(rnam))
                rownames(lambda.groups.mat)<-rnam
          } else {
            rnam<-lapply(nam,paste,sep="", collapse=",")
            if (nrow(lambda.groups.mat)==length(rnam))
                rownames(lambda.groups.mat)<-rnam
          }

      }

      if (exists("leg")) attr(lambda.groups.mat, which="objtable")<-leg

      ## worth matrix
      worth.groups.mat <- apply(lambda.groups.mat, 2, function(x) exp(2*x)/sum(exp(2*x)))
      attr(worth.groups.mat, which="objtable")<- attr(lambda.groups.mat, which="objtable")

      ## return
      switch(outmat,
         "lambda" = return(lambda.groups.mat),
         "worth" = return(worth.groups.mat),
         "est" = return(lambda.groups.mat),
         stop("     outmat must be either 'worth' or 'lambda'\n")
      )

## fitobj from patt*.fit
} else {

    envList<-obj$envList
    ncovpar<-envList$ncovpar

    Tmod<- regexpr("T",envList$resptype)>0               # check if time model

    if(!Tmod){                                           # if not a time model
         tpoints<-1
         nobj<-obj$envList$nobj
         npar<-(nobj - 1) * ncovpar
         lambda<-obj$coefficients[1:npar]
         lmat<-matrix(lambda,nrow=nobj-1)
         lmat<-rbind(lmat,rep(0,ncol(lmat)))
    } else {                                             # time model
         tpoints<-envList$tpoints
         nobj<-(envList$nitems-1)
         npar<-nobj*ncovpar*tpoints
         lambda<-obj$coefficients[1:npar]
         lmat<-matrix(lambda,nrow=nobj)
         lmat<-rbind(lmat,rep(0,ncol(lmat)))
         lmat<-matrix(lmat,ncol=tpoints)
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
    est<-matrix(group.est,nrow=nobj/tpoints)

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

} # end output from patt*fit

}
