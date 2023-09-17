runCircos <- function(ct_calls,
                      SV_file,
                      CNA_file,
                      outputdir,
                      genome.build="hg19",
                      circos_path="~/Applications/circos/current/bin/circos",
                      svclasses=c("deletion"="del",
                                  "tandemdup"="td",
                                  "head2headI"="h2h",
                                  "tail2tailI"="t2t",
                                  "translocation"="tra"))
{
    suppressWarnings(system(paste0("mkdir ",outputdir)))
    setwd(outputdir)
    ## ####################################################
    require(GenomicRanges)
    require(parallel)
    require(RColorBrewer)
    ## ####################################################
    svs <- read.table(SV_file,sep="\t",header=T)
    cnas <- read.table(CNA_file,sep="\t",header=T)
    samples <- unique(as.character(cnas[,1]))
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
    names(lSV) <- names(lCN) <- samples
    ## ####################################################
    classtocol <- function(class)
    {
        if(class==svclasses["deletion"])
            return("blue")
        if(class==svclasses["tandemdup"])
            return("red")
        if(class==svclasses["translocation"])
            return("grey")
        if(class==svclasses["head2headI"])
            return("dyellow")
        if(class==svclasses["tail2tailI"])
            return("lgreen")
        return("grey")
    }

    writeSegDup <- function(t,dirOut)
    {
        cols <- paste0("color=",sapply(t[,"svclass"],classtocol))
        nsvs <- cbind(paste0("hs",gsub("chr","",t[,"chr1"])),
                      t[,"pos1"],
                      t[,"pos1"]+1,
                      paste0("hs",gsub("chr","",t[,"chr2"])),
                      t[,"pos2"],
                      t[,"pos2"]+1,
                      cols)
        write.table(nsvs,
                    file=paste0(dirOut,"segdup.txt"),
                    quote=F,row.names=F,col.names=F)
    }

    reformatSF <- function(subclones, dirOut)
    {
        subclones$chr <- paste0("hs" ,gsub("chr","",subclones$chr))
        ## use allele specific copy number to offset the cicos plot highlights 1 CN = 50p difference (no max for now)
        subclones[is.na(subclones)] <- 0 ## not sure why needed but NA should not happen
        ntotfit <- (subclones$nMaj1+subclones$nMin1)## *subclones$frac1_A + (subclones$nMaj2_A+subclones$nMin2_A)*subclones$frac2_A -- the added fraction is for Battenberg calls but not all callers have this info
        nminfit <- (subclones$nMin1)##*subclones$frac1_A + (subclones$nMin2_A)*subclones$frac2_A
        ## format highlight file column and write output files
        circparams <- paste0("fill_color=orange,offset=",round(ntotfit*50),"p")
        df.subcl.tot <- data.frame(subclones[, c("chr", "start", "end")], circparams)
        circparams <- paste0("fill_color=black,offset=",round(nminfit*50),"p")
        df.subcl.min <- data.frame(subclones[, c("chr", "start", "end")], circparams)
        write.table(rbind(df.subcl.tot, df.subcl.min),
                    file=paste0(dirOut,"subcl.txt"),
                    sep = " ",
                    row.names = F,
                    col.names = F,
                    quote = F)
    }

    writeCTFlag <- function(ct_sample,dirOut)
    {
        starts <- ct_sample[,"start"]
        ends <- ct_sample[,"end"]
        chrs <- gsub("chr","",gsub("[pq]","",ct_sample[,"chrArm"]))
        tabflags <- cbind(paste0("hs",chrs),as.character(starts),as.character(ends))
        tabflags[tabflags=="hs23"] <- "hsX"
        tabflags[tabflags=="hs24"] <- "hsY"
        if(is.vector(tabflags))
            tabflags <- t(as.matrix(tabflags))
        write.table(rbind(tabflags[,],c(NULL,NULL,NULL)),
                    file=paste0(dirOut,"flaggedchromothripsis.txt"),
                    quote=F,row.names=F,col.names=F)
        return(NULL)
    }

    copyCircosFiles <- function(circosdir, dirOut)
    {
        system(paste0("cp ",circosdir,"mycircosDJ.conf ",dirOut))
        system(paste0("cp ",circosdir,"myideogram.position.conf ",dirOut))
        system(paste0("cp ",circosdir,"myideogram.conf ",dirOut))
        system(paste0("cp ",circosdir,"bands.conf ",dirOut))
        system(paste0("cp ",circosdir,"myideogram.label.conf ",dirOut))
        system(paste0("cp ",circosdir,"myticks.conf ",dirOut))
        system(paste0("cp ",circosdir,"karyotype.human.",genome.build,".txt ",dirOut))
    }

    cleanUpCircosFiles <- function(dirOut)
    {
        system(paste0("rm -f ",dirOut,"/mycircosDJ.conf "))
        system(paste0("rm -f ",dirOut,"/myideogram.position.conf "))
        system(paste0("rm -f ",dirOut,"/myideogram.conf "))
        system(paste0("rm -f ",dirOut,"/bands.conf "))
        system(paste0("rm -f ",dirOut,"/myideogram.label.conf "))
        system(paste0("rm -f ",dirOut,"/myticks.conf "))
        system(paste0("rm -f ",dirOut,"/karyotype.human.",genome.build,".txt "))
    }
    ## ####################################################
    circosdir <- paste0(system.file(package="CTCallR"),"/")
    ## ####################################################

    ## ####################################################
    samps <- unique(as.character(ct_calls[,"samplename"]))
    for(samp in samps)
    {
        print(paste0("Plotting sample: ",samp))
        subct <- ct_calls[ct_calls[,"samplename"]==samp,]
        ct_sample <- subct[subct[,"Calls"]=="Chromothripsis",]
        dirOut <- paste0(outputdir,"/",samp,"/")
        suppressWarnings(system(paste0("mkdir ",dirOut)))
        setwd(dirOut)
        copyCircosFiles(circosdir, dirOut)
        writeSegDup(lSV[[samp]],dirOut)
        reformatSF(lCN[[samp]],dirOut)
        cts <- writeCTFlag(ct_sample,dirOut)
        cmdCircos <- paste0(circos_path,
                            " -conf ",
                            dirOut,"/mycircosDJ.conf -outputdir ",
                            dirOut," -outputfile ",
                            samp,"_im.png > /dev/null 2>&1 ")
        try(system(cmdCircos,wait=T),silent=F)
        cleanUpCircosFiles(dirOut)
    }
    ## ####################################################
}
## ####################################################
