testorder <- function(svs)
{
    suppressWarnings(
        COR <- try(cor.test(svs[,"pos1"],
                        svs[,"pos2"],
                        met="sp")$p.value,silent=T))
    COR
}
