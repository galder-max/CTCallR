testclass <- function(svs, types=c("Deletion",
                                  "Head-to-head_inversion",
                                  "Tail-to-Tail_inversion",
                                  "Tandem_duplication"))
{
    require(EMT)
    classes <- svs[,"svclass"]
    observed <- as.vector(table(classes)[types])
    observed[is.na(observed)] <- 0
    suppressWarnings(Chi2 <- try(chisq.test(observed,p=rep(1/4,4))$p.value,silent=T))
    if(inherits(Chi2,"try-error")) Chi2 <- NA
    Chi2
}
