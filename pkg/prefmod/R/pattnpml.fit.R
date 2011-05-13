"pattnpml.fit" <-
 function(formula,  random=~1,
                    k = 1,
                    design,
                    tol = 0.5,
                    startp = NULL,
                    EMmaxit=500,               #
                    EMdev.change=0.001,        #
                    pr.it = FALSE
                    )
{
     call<-match.call()
     dat<-design
     pluginz<-NULL
     N<-nrow(dat)
     offset<-rep(0,N)
     weights<-rep(1,N)
     maxit <- EMmaxit       #
     conv <- EMdev.change  #
     prit <- pr.it
     if (!is.null(startp))
        if(!(length(startp) == k & all.equal(sum(startp),1)))
           stop("startp incorrectly specified!")

     RET <- alldistPC(formula=formula,
                    random = random,
                    family = poisson,
                    data=dat,
                    k = k,
                    random.distribution="np",
                    tol = tol,
                    offset,
                    weights,
                    pluginz=pluginz,
                    na.action,
                    EMmaxit=maxit,             #
                    EMdev.change=conv,         #
                    lambda=0,                  #
                    damp=FALSE,                #
                    damp.power=1,              #
                    spike.protect=0,           #
                    #sdev,
                    #shape,
                    plot.opt=0,
                    verbose=FALSE,
                    startp = startp,
                    pr.it = prit
     )
     RET$call <- call # rh 2010-24-04 obj$call should contain call to pattnpml.fit
     RET
}
