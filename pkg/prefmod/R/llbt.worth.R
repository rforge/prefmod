llbt.worth <- function(fitobj, outmat = "worth")
{

#### updated on 2010-12-21
#### can now handle fitobjects from llbtPCfit and
#### from gnm, provided llbt.design (also new version )was used
#### object specific covariates are possible

# check if either model results from gnm (using llbt.design)
# or from llbtPC.fit

llbtdesignmodel <- FALSE
llbtPCfitmodel <- ("llbtMod" %in% class(fitobj))
if (!llbtPCfitmodel){
   dterm <- deparse(fitobj$call$data)            ## unten nocheinmal
   x <- tryCatch(get(dterm),error = function(e) FALSE)
   if("llbtdes" %in% class(x)) llbtdesignmodel <- TRUE
}
if (llbtdesignmodel == llbtPCfitmodel)
     stop("Function must use result either from llbtPCfit or having used llbtdesign")


## fitobj from llbtdesign, gnm
if (llbtdesignmodel){

      ## subject covariates

      # extract subject covariates from eliminate
      eterm<-deparse(fitobj$call$eliminate)
      eterms<-unlist(strsplit(eterm,":"))
      if("mu" %in% eterms)
         eterms<-eterms[-which(eterms=="mu")]

      ### achtung laenge 0
      if(length(eterms>0)) {

          # extract subject covariates from formula
          fterm<-deparse(fitobj$call$formula)
          fterms<-unique(unlist(strsplit(eterm,"[ ()+*:~]")))
          subjcov.names<-intersect(eterms,fterms)
          ### achtung laenge 0

          # set up subj cov design
          levlist<-lapply(fitobj$model[subjcov.names],levels)
          maxlev<-sapply(levlist,length)
          subjdes<-prefmod:::gfac2(maxlev)
          colnames(subjdes)<-names(maxlev)

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

      # remove NA estimates from coefficients
      notna<-!is.na(coefs)
      if (length(notna)>0){
          sumMat<-sumMat[,notna,drop=F]
          coefs<-coefs[notna]
      }

      # remove unnecessary terms from estimates and from coefficients
      if (length(diffs)>0){
        dc<-which(names(coefs)==diffs)
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

## fitobj from llbtPC.fit
} else {

      nobj <- fitobj$envList$nobj

      # remove category and undecided parameters from lambda
      lambda <- fitobj$coefficients[fitobj$ofInterest]
      if(any(grep("^g[0-9]|g[0-9]$|^u$",names(lambda))))
          lambda <- lambda[-(grep("^g[0-9]|^u$", names(lambda)))]
      lambda <- ifelse(is.na(lambda),0,lambda)

      # initialise lambda matrix
      npar<-length(lambda)
      lambda.mat <- matrix(, nrow=nobj, ncol=npar/nobj)

      # row and colnames for output matrix
      nam.lambda <- names(lambda)
      rownames(lambda.mat) <- nam.lambda[1:nobj]
      nam.new<-nam.lambda
      for (i  in nam.lambda[1:nobj]) # remove obj names from nam.new
          # nam.new<-gsub(paste("^",i,"$|^",i,"[^0-9]:+|[:alnum:]*:?(",i,"[^0-9]{0,2})$",sep=""),"",nam.new)
          nam.new<-gsub(paste("^",i,"$|^",i,":",sep=""),"",nam.new)
      nn<-nam.new
      nl<-nam.lambda
      colnames(lambda.mat) <- unique(nam.new)


      # fill lambda matrix lambda.mat with correct entries

      # case with covs
      if(npar>nobj){

          # sort names for model terms according to columns of model matrix
          covmat<-unique(model.matrix(fitobj$envList$formel,data=fitobj$data))
          cnam<-colnames(covmat)
          nam.cov<-sapply(1:ncol(covmat), # sort variable names within interactions terms
              function(i){txt<-sort(unlist(strsplit(cnam[i],":"))); #
              paste(txt,collapse=":")}
          )

          A<-nam.lambda[1:nobj]
          B<-nam.cov
          B[1]<-""
          lo<-as.vector(t(outer(A,B,paste,sep=":")))
          lo<-sub(":$","",lo)
          # estimates in same order as in model matrix
          lambda<-lambda[match(lo,names(lambda))]

          lambda.mat<-matrix(lambda,nr=nobj,b=T)

          colnames(lambda.mat) <- unique(nam.new)
          rownames(lambda.mat) <- nam.lambda[1:nobj]

          ## calculate sums of estimates accordig to covariate groups
          lambda.vec <- as.vector(lambda.mat)
          lambda.groups.mat <- matrix((covmat %x% diag(nobj)) %*% lambda.vec, nrow=nobj)


          ## labels for cov groups
          x<-fitobj$envList$model.covs
          xx<-mapply(function(x,y)paste(x,y,sep=""), colnames(x),data.frame(x))
          gr.labels <-apply(xx,1,paste,collapse=":")
          colnames(lambda.groups.mat) <- gr.labels


          ## obsolete since option obj.names removed 2011-01-03
          #if (is.null(obj.names))
          #   rownames(lambda.groups.mat) <- nam.lambda[1:nobj]
          #else
          #   rownames(lambda.groups.mat) <- obj.names
          rownames(lambda.groups.mat) <- nam.lambda[1:nobj]

          ## worth matrix
          worth.groups.mat <- apply(lambda.groups.mat, 2, function(x) exp(2*x)/sum(exp(2*x)))

      # case w/o covs
      } else {
      lambda.mat <- matrix(lambda, nrow = nobj, dimnames = list(names(lambda)[1:nobj], "estimate"))
      lambda.groups.mat <- lambda.mat
      worth.groups.mat <- apply(lambda.groups.mat, 2, function(x) exp(2*x)/sum(exp(2*x)))
      }
}


## return
switch(outmat,
   "lambda" = return(lambda.groups.mat),
   "worth" = return(worth.groups.mat),
   "est" = return(lambda.groups.mat),
   stop("     outmat must be either 'worth' or 'lambda'\n")
)

}
