doPCF <-
function(nt, minEvents)
{
    require(copynumber)
    diffs <- c(nt[1,"end"]-nt[1,"start"],diff(nt[,"end"]))
    res <- pcf(data.frame(chr=rep(nt$chr[1],length(diffs)),
                          positions=nt[,"end"],sample1=diffs),
               verbose=F,
               kmin=minEvents)
    res
}
