treatPCF <-
function(pcf,
                     flagDensity=150000000/50,
                     minEvents=5)
{
    npcf <- pcf
    whichFlaggedP <- which(pcf$mean<flagDensity & pcf$n.probes>=minEvents & pcf$arm=="p")
    if(length(whichFlaggedP)>0)
    {
        whichFlaggedP <- min(whichFlaggedP):max(whichFlaggedP)
        if(length(whichFlaggedP)>1)
            npcf <- npcf[-c(whichFlaggedP[2:length(whichFlaggedP)]),]
        i <- whichFlaggedP[1]
        npcf[i,"start.pos"] <- pcf[min(whichFlaggedP),"start.pos"]
        npcf[i,"end.pos"] <- pcf[max(whichFlaggedP),"end.pos"]
        npcf[i,"n.probes"] <- sum(pcf[whichFlaggedP,"n.probes"])
        wfp <- whichFlaggedP
        npcf[i,"mean"] <- sum(pcf[wfp,"mean"]/1000*pcf[wfp,"n.probes"])/sum(pcf[wfp,"n.probes"])*1000
    }
    pcf <- npcf
    whichFlaggedP <- which(pcf$mean<flagDensity & pcf$n.probes>=minEvents & pcf$arm=="q")
    if(length(whichFlaggedP)>0)
    {
        whichFlaggedP <- min(whichFlaggedP):max(whichFlaggedP)
        if(length(whichFlaggedP)>1)
            npcf <- npcf[-c(whichFlaggedP[2:length(whichFlaggedP)]),]
        i <- whichFlaggedP[1]
        npcf[i,"start.pos"] <- pcf[min(whichFlaggedP),"start.pos"]
        npcf[i,"end.pos"] <- pcf[max(whichFlaggedP),"end.pos"]
        npcf[i,"n.probes"] <- sum(pcf[whichFlaggedP,"n.probes"])
        wfp <- whichFlaggedP
        npcf[i,"mean"] <- sum(pcf[wfp,"mean"]/1000*pcf[wfp,"n.probes"])/sum(pcf[wfp,"n.probes"])*1000
    }
    return(npcf)
}
