deriveTable <-
function(allCPs)
{
    tt <- NULL
    for(i in 1:length(allCPs))
    {
        cat(".")
        for(j in 1:length(allCPs[[i]]))
        {
            cat("+")
            if(allCPs[[i]][[j]]$flagged & !is.na(allCPs[[i]][[j]]$flagged))
            {
                if(allCPs[[i]][[j]]$flaggedDensity)
                {
                    vv <- tablise(allCPs[[i]][[j]],j,names(allCPs)[i])
                    tt <- rbind(tt,vv)
                }
            }
        }
    }
    colnames(tt) <- c("samplename","start","end",
                      "chrArm",
                      "density",
                      "NbreakpointsCNA",
                      "p.KS_exp",
                      "CNA.states",
                      "CNA.sizes",
                      "coverage.mode",
                      "coverage.max")
    tt
}
