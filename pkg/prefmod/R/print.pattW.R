print.pattW<-function(x,...)
{
    x<-unclass(x)
    cat("\nWorthmatrix:\n")
    if(ncol(x)>1) cat("\n")  # blank line only if there are groups
    print(round(x,digits=4))
    cat("\n")
}
