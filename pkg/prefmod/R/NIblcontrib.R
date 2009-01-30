# calculates likelihood contributions for MNAR
#
NIblcontrib<-function(obj,lambda,X,nobj,ENV)
{
      naidx<-1-as.numeric(obj$notnaidx)  # R=1 missing, R=0 observed

      R<-matrix(rep(naidx,2^ENV$ncomp),nr=2^ENV$ncomp,byrow=T)
      RBstar<-R %*%abs(pcdesign(nobj))  # alpha_i + alpha_j
      #RBstar<-R %*%(pcdesign(4))    # alpha_i - alpha_j

      YRBstar<-do.call(cbind,lapply(1:(nobj),function(i) RBstar[,i]*ENV$Y[,i])) # betas

      #####    nonresponse model alpha_i+alpha_j/alpha_i-alpha_j
      XX<-X                              # only lambdas
 ##     if (ENV$MISalpha) XX<-cbind(X,RBstar[,1:nobj])          # lambdas, alpha_i
      if (ENV$MISalpha) XX<-cbind(X,RBstar[,ENV$Malph])          # lambdas, alpha_i
 ##     if (ENV$MISbeta)  XX<-cbind(X,RBstar[,1:nobj],YRBstar)  # lambdas, alpha_i. beta_i
      if (ENV$MISbeta)  XX<-cbind(X,RBstar[,ENV$Malph],YRBstar[,ENV$Mbeta])  # lambdas, alpha_i. beta_i


      #####    nonresponse model common alpha
      if (ENV$MIScommon) {
         r<-sum(naidx)
         XX<-cbind(X,r)                                       # lambdas, alpha+alpha
      }


      # add dependencies if ia=TRUE
      if(ENV$ia) XX<-cbind(XX,ENV$XI)

      patt.all <- exp(XX %*% lambda)
      patt<-tapply(patt.all,obj$s,sum)
      patt<-patt/ENV$nrm.all  #  divide by normalizing constant
      ll<-sum(obj$counts*log(patt))


      # log likelihood for saturated model
      p.cnts<-obj$counts/sum(obj$counts)
      fl<-sum(log(p.cnts[p.cnts>0])*obj$counts[obj$counts>0])

      ###RET<-list(ll1=ll1,fl=fl,summ.patt.all=summ.patt.all,summ.cnts=summ.cnts)
      RET<-list(ll=ll,fl=fl)
      RET
}
