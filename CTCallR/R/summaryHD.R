summaryHD <-
function(cn)
{
    hd <- (cn$nMaj1==0 & cn$nMin1==0) | (cn$nMaj2==0 & cn$nMin2==0)
    sizes <- sum(cn[hd,"end"]/1000000-cn[hd,"start"]/1000000,na.rm=T)
}
