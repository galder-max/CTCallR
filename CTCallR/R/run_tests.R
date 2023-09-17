run_tests <- function(ct, lSV, samples, svtypes)
{
    chr <- gsub("chr","",gsub("p","",gsub("q","",ct[4])))
    if(chr=="23") chr <- "X"
    svs <- lSV[[which(samples==ct["samplename"])]]
    require(GenomicRanges)
    grSV1 <- GRanges(svs[,"chr1"],IRanges(svs[,"pos1"],svs[,"pos1"]))
    grSV2 <- GRanges(svs[,"chr2"],IRanges(svs[,"pos2"],svs[,"pos2"]))
    grCT <- GRanges(chr,
                    IRanges(as.numeric(ct[2]),as.numeric(ct[3])))
    ov1 <- findOverlaps(grSV1,grCT)
    ov2 <- findOverlaps(grSV1,grCT)
    keep <- unique(c(queryHits(ov1),queryHits(ov2)))
    svs <- svs[keep,]
    list(testclass(svs, svtypes),
         testorder(svs[svs[,"chr1"]==chr & svs[,"chr2"]==chr,]),
         length(keep))
}
