returnCov <-
function(nb_cna)
{
    b <- (.5-100/50*.8)/(1-100/50)
    a <- (.8-b)/50
    ## linear cap depending on number fragments; max at 1
    cov <- a*nb_cna+b
    cov[cov>1] <- 1
    cov
}
