mytry <-
function(x,retVal=NA,...)
{
    kk <- try(x,silent=T,...)
    if(inherits(kk,"try-error")) return(retVal)
    kk
}
