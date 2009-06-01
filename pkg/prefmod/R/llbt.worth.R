llbt.worth <- function(obj, obj.names=NULL, outmat="worth"){

   #if(!(class(obj)[1] %in% c("glm","gnm"))) stop
   if(!("llbtMod" %in% class(obj)))
        stop("function only for objects of class llbtMod (see help for llbtPC.fit)")

   nobj <- obj$envList$nobj

   mtype <- ifelse(is.null(obj$ofInterest),"glm","gnm")


   if(mtype == "gnm"){
       lambda <- obj$coefficients[obj$ofInterest]
       if(any(grep("^g[0-9]|g[0-9]$|^u$",names(lambda))))
           lambda <- lambda[-(grep("^g[0-9]|^u$", names(lambda)))]
       lambda <- ifelse(is.na(lambda),0,lambda)

       npar<-length(lambda)
       lambda.mat <- matrix(, nrow=nobj, ncol=npar/nobj)

       # row and colnames for output matrix
       nam.lambda <- names(lambda)
       rownames(lambda.mat) <- nam.lambda[1:nobj]
       nam.new<-nam.lambda
       for (i  in nam.lambda[1:nobj])
#           nam.new<-gsub(paste("^",i,"$|^",i,"[^0-9]:+|[:alnum:]*:?(",i,"[^0-9]{0,2})$",sep=""),"",nam.new)
           nam.new<-gsub(paste("^",i,"$|^",i,":",sep=""),"",nam.new)
       nn<-nam.new
       nl<-nam.lambda
       colnames(lambda.mat) <- unique(nam.new)

       # fill lambda matrix with correct entries
       npar<-length(lambda)

       # case with covs
       if(npar>nobj){
           lambda.mat <- matrix(, nrow=nobj, ncol=npar/nobj)

           # reference group into first column
           lambda.mat[,1] <- lambda[1:nobj]

           # others according to their labels
           lambda.rest<-lambda[-(1:nobj)]
           nn2<-nn[-(1:nobj)]
           nn2<-nl[-(1:nobj)]
           lambda.mat.rest<-matrix(lambda.rest,b=T,nrow=nobj)
           lambda.mat[,2:ncol(lambda.mat)]<-lambda.mat.rest
       colnames(lambda.mat) <- unique(nam.new)
       rownames(lambda.mat) <- nam.lambda[1:nobj]

       ## calculate sums accordig to covariates
       lambda.vec <- as.vector(lambda.mat)
       covmat <- unique(obj$envList$covdesmat)

       lambda.groups.mat <- matrix((covmat %x% diag(nobj)) %*% lambda.vec, nrow=nobj)


       ## labels for cov groups
       x<-obj$envList$model.covs
       xx<-mapply(function(x,y)paste(x,y,sep=""), colnames(x),data.frame(x))
       gr.labels <-apply(xx,1,paste,collapse=":")
       colnames(lambda.groups.mat) <- gr.labels

       if (is.null(obj.names))
          rownames(lambda.groups.mat) <- nam.lambda[1:nobj]
       else
          rownames(lambda.groups.mat) <- obj.names

       ## worth matrix
       worth.groups.mat <- apply(lambda.groups.mat, 2, function(x) exp(2*x)/sum(exp(2*x)))

       # case w/o covs
       } else {
       lambda.mat <- matrix(lambda, nrow = nobj, dimnames = list(names(lambda)[1:nobj], "estimate"))
       lambda.groups.mat <- lambda.mat
       worth.groups.mat <- apply(lambda.groups.mat, 2, function(x) exp(2*x)/sum(exp(2*x)))
       }

   } else if(mtype == "glm"){
   } else {
   }

   switch(outmat,
      "lambda" = return(lambda.groups.mat),
      "worth" = return(worth.groups.mat),
      "est" = return(lambda.mat),
      stop("     outmat must be either 'worth' or 'lambda'\n")
   )
   #list(lam=lambda.mat, lamsum=lambda.groups.mat, w=worth.groups.mat)
}
