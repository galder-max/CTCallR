# CTCallR

![](
"CTCallR - A stringent ChromoThripsis Caller from human WGS in R")

The CTCallR is an R package destined at calling chromothripsis from
whole-genome sequencing data. It applies a series of tests proposed in
(Korbel & Campbell, Cell,
2013)[https://doi.org/10.1016/j.cell.2013.02.023] to annotate
high-density-breakpoint regions as "clusters of breakpoints" or
(likely) "chromothripsis".

It expects the user to have run a copy-number and structural variant
caller and format their output properly.


## Install

You can install the package with devtools in an R session:

> devtools::install_github("VanLoo-lab/CTCallR", build_opts = c("--no-build-vignettes"))

## Dependencies 

The package relies on a few dependencies, nameley the R packages: vcf,
MASS, copynumer, parallel, EMT.

It also relies on the (circos)[http://circos.ca/] software to plot the results.


## Other CT caller

Please see also (ShatterSeek)[https://github.com/parklab/ShatterSeek], for another chromothripsis caller used in The
Pan-Cancer Analysis of Whole Genomes (PCAWG) study. 


## Example run

To test CTCallR, we propose example data from 3 sarcomas of the ICGC PCAWG
UK cohort.


## Usage

See documentation for example run. In an R session, type:

> ?CTCallR


