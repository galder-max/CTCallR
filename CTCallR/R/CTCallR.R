CTCallR <- function(CNA_file,
                    SV_file,
                    outputdir,
                    plot_results=TRUE,
                    outputfile="./chromothripsis_annotated_calls.txt",
                    svclasses=c("deletion"="del",
                                "tandemdup"="td",
                                "head2headI"="h2h",
                                "tail2tailI"="t2t",
                                "translocation"="tra"),
                    prandomclass_threshold=0.01,##absence of rejection so we are stricter by default
                    psegmentsize_threshold=0.05, ##needs a rejection so we put the p-value high not to be too strict
                    min_cna_breakpoints=15, ##min number of CNA breakpoints for high density
                    min_sv_breakpoints=15, ##min number of SV breakpoints for high density
                    flagDensity=150000000/50, ##min local density of breakpoints to flag region for high density
                    minEvents=5, ##minimum number of events for segmentation
                    genome.build="hg19",
                    circos_path="~/Applications/circos/current/bin/circos",
                    MC.CORES=5)
{
    require(vcd)
    require(MASS)
    require(copynumber)
    require(parallel)
    require(EMT)
    ## ############################################
    N <- NULL
    svtypes <- svclasses[!names(svclasses)=="translocation"]
    ## ############################################
    ## ############################################
    ## Load SVs and CNAs
    svs <- read.table(SV_file,sep="\t",header=T)
    cnas <- read.table(CNA_file,sep="\t",header=T)
    ## ############################################
    ## ############################################
    ## Unique samples
    samples <- unique(as.character(cnas[,1]))
    ## ############################################
    ## list of CNA and SV, columns renamed
    lCN <- lapply(samples,function(x)
    {
        tmp <- cnas[as.character(cnas[,1])==x,]
        colnames(tmp) <- c("samplename","chr","start","end",
                           "nMaj1","nMin1","frac1")
        tmp[order(tmp[,2],tmp[,3],decreasing=F),]
    })
    lSV <- lapply(samples,function(x)
    {
        tmp <- svs[as.character(svs[,1])==x,]
        colnames(tmp) <- c("samplename",
                           "chr1","pos1",
                           "chr2","pos2",
                           "svclass")
        tmp
    })
    ## ############################################
    ## ############################################
    ## quick sanity check on size of 0+0 regions
    hd <- sapply(lCN,summaryHD)
    print("Homozygous deletions:")
    print(summary(hd))
    ## ############################################
    ## ############################################
    ## list of flagged high density breakpoint regions
    ## also annotates p-value for test for distribution of segment sizes
    print("Process to flagging regions sample by sample")
    allCPs <- mclapply(lCN,function(t)
    {
        cat(".")
        CPs <- lapply(c(1:22,"X"),function(x) {
            CP(t[as.character(t$chr)==x,], flagDensity, minEvents)
        })
    },mc.cores=MC.CORES)
    ## ############################################
    names(allCPs) <- samples
    ## ############################################
    ## make table from flagged regions
    tt <- deriveTable(allCPs)
    ## ############################################
    ## ############################################
    ## test regions for random SV classes/types; and random mate position order (not used)
    print("Run statistical tests on CNA and SVs across samples")
    alltest <- apply(tt,1,function(x) run_tests(x, lSV, samples, svtypes)) ## all tests are performed: expected warnings for Chi^2 approx.
    pClass <- sapply(alltest,function(x) x[[1]]) ## p-values for random classes
    pOrder <- sapply(alltest,function(x) x[[2]]) ## p-values for random mate order
    nbSVs <- sapply(alltest,function(x) x[[3]]) ## number of breakpoints
    ## ############################################
    ## Augment table with p-values; expected min. region size covered by at most 3 copy number states;
    ## and number of breakpoints;
    print("Annotate flagged regions")
    tt <- cbind(tt,pClass,pOrder)
    colnames(tt)[(ncol(tt)-1):ncol(tt)] <- c("pRandomClass","pRandomOrder")
    tt <- cbind(tt,returnCov(as.numeric(tt[,"NbreakpointsCNA"])))
    colnames(tt)[ncol(tt)] <- c("expectedCov")
    tt <- cbind(tt,nbSVs)
    colnames(tt)[ncol(tt)] <- c("nbSVs")
    ## ############################################
    ## ############################################
    ## Goes through flagged regions and annotations and annotates as "Chromothripsis" or "Cluster"
    print("Call chromothripsis from all flags and tests")
    isChromothripsis <- apply(tt,1,function(x)
    {
        mode <- as.numeric(x["coverage.mode"])>as.numeric(x["expectedCov"])
        max <- as.numeric(x["coverage.max"])>as.numeric(x["expectedCov"])
        class <- as.numeric(x["pRandomClass"])>prandomclass_threshold
        expt <- as.numeric(x["p.KS_exp"])<psegmentsize_threshold
        nbsv <- as.numeric(x["NbreakpointsCNA"])>min_cna_breakpoints
        nbcna <- as.numeric(x["nbSVs"])>min_sv_breakpoints
        (mode|max)&class&expt&(nbsv|nbcna)
    })
    tt <- cbind(tt,ifelse(isChromothripsis,"Chromothripsis","Cluster"))
    colnames(tt)[ncol(tt)] <- "Calls"
    ## ############################################
    ## ############################################
    ## Make data.frame (not sure why I put this there)
    print("Put it all together in a table format")
    tt <- as.data.frame(tt)
    tt$start <- as.numeric(as.character(tt$start))
    tt$end <- as.numeric(as.character(tt$end))
    tt$p.KS_exp <- as.numeric(as.character(tt$p.KS_exp))
    tt$pRandomClass <- as.numeric(as.character(tt$pRandomClass))
    suppressWarnings(tt$pRandomOrder <- as.numeric(as.character(tt$pRandomOrder)))
    rownames(tt) <- NULL
    ## ############################################
    ## save table of annotated regions and calls
    print("Write to disk")
    write.table(tt,
                file=outputfile,
                sep="\t",
                col.names=T,
                row.names=F,
                quote=F)
    ## ############################################
    if(plot_results)
    {
        print("Plot results")
        runCircos(ct_calls=tt,
                  SV_file=SV_file,
                  CNA_file=CNA_file,
                  outputdir=outputdir,
                  genome.build=genome.build,
                  circos_path=circos_path,
                  svclasses=svclasses)
    }
    ## ############################################
    print("Return table")
    return(tt)
    ## ############################################
}
