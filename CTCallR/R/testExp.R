testExp <- function(chrCN,N=NULL)
{
    require(vcd)
    require(MASS)
    sizes <- diff(chrCN[,"end"])
    if(length(sizes)<4) stop("number of breakpoints too low")
    if(any(sizes<0))
    {
        print(chrCN)
        stop("non-ordered chromosomal positions")
    }
    fit <- fitdistr(sizes,"exponential")
    ks.test(sizes,pexp,rate=fit$estimate)$p.value
}
