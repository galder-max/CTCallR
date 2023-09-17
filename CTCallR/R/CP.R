CP <-
function(chrCN,
               flagDensity,
               minEvents)
{
    if(nrow(chrCN)>minEvents)
    {
        pcf <- doPCF(chrCN,minEvents=minEvents)
        pcf <- treatPCF(pcf,flagDensity,minEvents)
        starts <- cumsum(c(1,pcf$n.probes[-c(nrow(pcf))]))
        ends <- cumsum(pcf$n.probes)
        whichFlagged <- which(pcf$mean<flagDensity & pcf$n.probes>=minEvents)
        tE <- mytry(testExp(chrCN,N=NULL))
        tERegion <- lapply(whichFlagged,function(x)
        {
            res <- try(testExp(chrCN[starts[x]:ends[x],],N=N),silent=T)
            if(inherits(res,"try-error")) return(NULL)
            res
        })
        flaggedtE <- tE<0.05
        flaggedDensity <- length(whichFlagged)>0
        flagged <- flaggedtE | flaggedDensity
        copyNstates <- getCNstates(pcf,whichFlagged,chrCN)
        return(list(starts=chrCN[starts[whichFlagged],"start"],
                    ends=chrCN[ends[whichFlagged],"end"],
                    arms=pcf[whichFlagged,"arm"],
                    densities=pcf$mean[whichFlagged],
                    Nevents=pcf$n.probes[whichFlagged],
                    cnStates=copyNstates,
                    testExpChromosome=tE,
                    testExpRegion=tERegion,
                    flaggedtE=flaggedtE,
                    flaggedDensity=flaggedDensity,
                    flagged=flagged))
    }
    else
        return(list(flagged=FALSE,
                    minEventsReached=FALSE))
}
