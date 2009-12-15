# calculates normalizing constant for all Y,R patterns

normalize<-function(lambda,X,nobj,ENV){

   nrm<-function(r,eta,Bstar,alpha,betas,Xundec,XIthet,ENV){   # 7.12.09  Xundec added
      if(ENV$MISalpha){
        rBstar<-r %*% Bstar                  # one row of RBstar, identical for all rows within block
        eta<-eta + as.vector(rBstar[,ENV$Malph]%*%alpha) # scalar, based on one row of RBstar,
                                             # identical for all rows within block
      }
      if (ENV$MISbeta){
        YrBstar<-do.call(cbind,lapply(1:nobj,function(x)ENV$Y[,x]*rBstar[x])) # design matrix for betas for current block
        eta<-eta+YrBstar[,ENV$Mbeta]%*%betas
      }

      if(ENV$undec){                 # 7.12.09
        eta<-eta+Xundec
      }
      if(ENV$ia){
        eta<-eta+XIthet
      }

      if(ENV$MIScommon){
        eta<-eta+as.vector(r)*alpha
      }

      nrm.bl<-sum(exp(eta))
      nrm.bl
   }# end nrm()


   #alpha<-betas<-theta<-NULL
   # lam(bda)s, alphas, betas from parametervector lambda
   lam<-lambda[1:(nobj-1)]
   if(ENV$MISalpha)
     alpha<-lambda[ENV$paridx==2]
   if(ENV$MISbeta)
     betas<-lambda[ENV$paridx==3]
   if(ENV$MIScommon){
     alpha<-lambda[ENV$paridx==4]
   }
   if(ENV$undec){                       # 7.12.09
     undecU<-lambda[ENV$paridx==6]
     Xundec<-ENV$U * undecU
   }
   if(ENV$ia){
     theta<-lambda[ENV$paridx==5]
     XIthet<-ENV$XI%*%theta
   }

   eta<-X%*%lam # YBlambda

   R<-patternmat2(nobj)
   if(ENV$MIScommon)
     R<-matrix(rowSums(R),nc=1)

   Bstar<-abs(pcdesign(nobj))  # alpha_i + alpha_j
   nrm.bl.vec<-apply(R,1,nrm,eta,Bstar,alpha,betas,Xundec,XIthet,ENV)        # 7.12.09  Xundec added
   nrm.all<-sum(nrm.bl.vec)
   nrm.all
}
