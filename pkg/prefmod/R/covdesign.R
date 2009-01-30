# design matrix for subject covs
covdesign<-function(formel,covs)
{
   # number of factor levels for terms from formula
   vars<-attr(terms(formel),"term.labels")
   order<-attr(terms(formel),"order")
   vars<-vars[order==1]
   levs<-apply(as.matrix(covs[,vars]),2,max) # number of factor levels


   maineffects<-as.data.frame(gfac2(levs)) #
   maineffects<-data.frame(apply(maineffects,2,factor))
   names(maineffects)<-vars

   form<-formula(formel)
   mmat<-model.matrix(form,data=maineffects)  # intercept included
   as.matrix(unique(mmat))
}
