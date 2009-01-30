# split data according to subject covariates
splitCovs<-function(dat,covs,formel,elim,ENV)
{
   if(is.null(covs)){
      formel<-~1                 # remove formulae if no covariates
      elim<-~1
   }

   if(elim=="~1") elim<-formel

   if(elim=="~1"){               # subject covariates
      covdesmat<-matrix(1,1)
      elimdesmat<-covdesmat
      covs<-as.matrix(rep(1,nrow(dat)),nc=1)
   } else {
      if (formel=="~1"){
           covdesmat<-matrix(1,1)
           colnames(covdesmat)<-""
      } else {
           covdesmat<-covdesign(formel,covs)
      }
      elimdesmat<-covdesign(elim,covs)
   }

   all.term.labels<-attr(terms(elim),"term.labels")
   order<-attr(terms(elim),"order")
   all.frml.term.labels<-attr(terms(formel),"term.labels")
   order.frml<-attr(terms(formel),"order")

   ENV$covlevels<-NULL

   if (length(order)>0){

       maineffect.terms<-all.term.labels[order==1]
       maineffect.terms.frml<-all.frml.term.labels[order.frml==1]
       if(length(maineffect.terms.frml)>0){                                        # only if cov terms in formel
         ENV$covlevels<-apply(as.matrix(covs[,maineffect.terms.frml]),2,max)
         names(ENV$covlevels)<-maineffect.terms.frml
       }

       # split data according to subj cov groups in eliminate formula
       cList<-split(dat,covs[,maineffect.terms])

       # remove terms in elim design matrix not in formel design matrix
       enam<-colnames(elimdesmat)                                                # had to be introduced
       cnam<-colnames(covdesmat)                                                 # since terms in
       nam.elim<-sapply(1:ncol(elimdesmat),                                      # interactions
                         function(i){txt<-sort(unlist(strsplit(enam[i],":")));   # may have different
                                     paste(txt,collapse=":")}                    # order in elim and formel
                       )                                                         # correct in R 2.6.x and
       nam.cov<-sapply(1:ncol(covdesmat),                                        # maybe in 2.7.0
                         function(i){txt<-sort(unlist(strsplit(cnam[i],":")));   #
                                     paste(txt,collapse=":")}                    #
                       )                                                         #
                                                                                 #
       model.terms<-nam.elim %in% nam.cov                                        #
       #model.terms<-colnames(elimdesmat) %in% colnames(covdesmat)          replaced

       colnames(elimdesmat)<-nam.elim                                            #

       for (i in 1:length(cList))
           cList[[i]]<-c(cList[i],list(cov=elimdesmat[i,model.terms]))
       ENV$covdesmat  <-elimdesmat[,model.terms]

   } else {

       cList<-split(dat,as.factor(rep(1,nrow(dat))))
       cList[[1]]<-c(cList[1],list(cov=1))
       ENV$covdesmat<-covdesmat
   }

   ENV$ncovpar<-ncol(covdesmat)
   ENV$formel<-formel
   ENV$elim<-elim
   cList

}
