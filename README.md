# CTCallR


The CTCallR is an R package destined at calling chromothripsis from
whole-genome sequencing data. It applies a series of tests proposed in
[Korbel & Campbell, Cell,
2013](https://doi.org/10.1016/j.cell.2013.02.023) to annotate
high-density-breakpoint regions as "clusters of breakpoints" or
(likely) "chromothripsis".

It expects the user to have run a copy-number and structural variant
caller and format their output properly.

An example output figure is shown below:

![](https://github.com/galder-max/CTCallR/blob/main/CTCallR/f393bb05-c737-4cc3-e040-11ac0d48452a_im.svg
"CTCallR - A stringent ChromoThripsis Caller from human WGS in R")


## Install

You can install the package with devtools in an R session:

> devtools::install_github("VanLoo-lab/CTCallR", build_opts = c("--no-build-vignettes"))

## Dependencies 

The package relies on a few dependencies, nameley the R packages: vcf,
MASS, copynumer, parallel, EMT.

It also relies on the [circos](http://circos.ca/) software to plot the results.


## Other CT caller

Please see also [ShatterSeek](https://github.com/parklab/ShatterSeek), for another chromothripsis caller used in The
Pan-Cancer Analysis of Whole Genomes (PCAWG) study. 


## Required input and their format
<br>
The main functions require two input files: copy-number calls across
samples, as well as structural variant calls across samples. 
<br>
Ideally, copy-number calls should have been
informed with the structural variants to identify the copy-number
breakpoints, as is done in Battenberg.
<br>
Example files are given in the example_input directory for both input
file types.<br>
<br>
For the copy number calls a tsv file, where the rows represent genomic
segments, with the following columns in
that order:<br>
<br>
samplename: the name of the sample in which this segment has been characterized<br>
chr: the chromosome of that segment<br>
startpos	: the start position of that segment<br>
endpos	: the end position of that segment<br>
nMaj1_A	: the number of copies of the major allele<br>
nMin1_A	: the number of copies of the minor allele<br>
frac1_A: the fraction subclonal (leave as 1 if no subclonality
information)<br>
<br>
For the structural variant calls, a tsv file, where the rows represent
pairs of structural variant breakpoints, with columns in that order:<br>
<br>
samplename : the name of the sample in which this sv has been characterized<br>
chr1	: the first chromosome of the pair<br>
pos1	: the position on that first chromosome of the pair<br>
chr2	: the second chromosome of the pair<br>
pos2	: the position on that second chromosome of the pair<br>
class   : the type of sv (can be extracted from most sv callers with
strand and mapping info), either of the type head-to-head inversion,
tail-to-tail inversion, tandem duplication, deletion or translocation.<br>




## Example run

To test CTCallR, we propose example data from 3 sarcomas of the ICGC PCAWG
UK cohort.


## Usage

See documentation for example run. In an R session, type:

> ?CTCallR


