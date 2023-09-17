detach("package:CTCallR",unload=T)
library(CTCallR)

CNA_file <- "CTCallR/CTCallR/example_input/allCNAs.txt"

SV_file <- "CTCallR/CTCallR/example_input/allSVs.txt"

setwd("/Users/tarabim/Desktop/Projects/Chromothripsis/github/")

if(F) debug(CTCallR)

ct <- CTCallR(CNA_file,
              SV_file,
              outputfile = "./chromothripsis_annotated_calls.txt",
              outputdir="/Users/tarabim/Desktop/Projects/Chromothripsis/github",
              circos_path="~/Desktop/Utils/circos-0.69-9/bin/circos")
