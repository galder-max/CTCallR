getCNstates <-
function(pcf,whichFlagged,chrCN,maxStates=3)
{
    starts <- c(0,cumsum(pcf$n.probes[-c(nrow(pcf))]))+1
    ends <- c(cumsum(pcf$n.probes))
    nStates <- list()
    for(i in whichFlagged)
    {
        wF <- as.logical(1:nrow(chrCN)*0)
        wF[starts[i]:ends[i]] <- TRUE
        states <- paste(chrCN[wF,"nMaj1"],chrCN[wF,"nMin1"],sep="-")
        sizes <- chrCN[wF,"end"]-chrCN[wF,"start"]
        isSubclonal <- chrCN[wF,"frac1"]
        ssizes <- tapply(sizes,as.factor(states),sum)/1000
        nsegs <- tapply(sizes,as.factor(states),length)
        fsizes <- ssizes/sum(ssizes)
        nS <- length(ssizes)
        if(nS<maxStates) maxStates <- nS
        nStates[[i]] <- list(states=states,
                             sizes=sizes,ssizes=ssizes,
                             isSubclonal=isSubclonal,
                             subclonalFrac=sum(sizes/100*as.numeric(isSubclonal!=1),na.rm=T)/sum(sizes/100,na.rm=T),
                             fsizes=fsizes,
                             nsegs=nsegs,
                             coveragemaxStates.modeSegs=sum(fsizes[order(nsegs,decreasing=T)][1:maxStates]),
                             coveragemaxStates.maxCov=sum(sort(fsizes,decreasing=T)[1:maxStates]))
    }
    nStates
}
