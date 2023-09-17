\name{CTCallR}
\alias{CTCallR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
CTCallR is the main function of the chromothripsis caller. 
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function loads a copy-number profile and a structural-variant call
file from disk, then goes through all segments and genomic regions to
identify regions with high-density of breakpoints and test their
statistical significance for different chromothripsis hallmarks. It will
generate a final table with all the regions and their annotation.
}
\usage{
CTCallR(CNA_file, SV_file, outputfile = "./chromothripsis_annotated_calls.txt", prandomclass_threshold = 0.01, psegmentsize_threshold = 0.05, min_cna_breakpoints = 15, min_sv_breakpoints = 15, flagDensity = 1.5e+08/50, minEvents = 5, MC.CORES = 5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{CNA_file}{
%%     ~~Describe \code{CNA_file} here~~
character string. path to the copy-number profile. this must be
formatted in a Battenberg-like fashion.
}
  \item{SV_file}{
%%     ~~Describe \code{SV_file} here~~
character string. path to the structural variant calls. this must be
formatted similarly to pcawg sv calls and must include the svclass for
each sv (head-to-head inv, tail-to-tail inv, deletion or tandem duplication).
}
  \item{outputfile}{
%%     ~~Describe \code{outputfile} here~~
character string. full path of the file where to save the final table with the regions and their annotation.
}
  \item{prandomclass_threshold}{
%%     ~~Describe \code{prandomclass_threshold} here~~
numeric value. p-value threshold for the test of uniformity of SV
classes. expects absence of rejection for CT.
}
  \item{psegmentsize_threshold}{
    %%     ~~Describe \code{psegmentsize_threshold} here~~
    numeric value. p-value threshold for the test of exponential distribution
    of segment sizes, which is expected under a sequential breakpoint
    model. expects rejection for CT.
}
  \item{min_cna_breakpoints}{
    %%     ~~Describe \code{min_cna_breakpoints} here~~
    numeric value. minimum number of copy-number breakpoints to consider CT.
    
}
  \item{min_sv_breakpoints}{
    %%     ~~Describe \code{min_sv_breakpoints} here~~
    numeric value. minimum number of SV breakpoints to consider CT.

}
  \item{flagDensity}{
    %%     ~~Describe \code{flagDensity} here~~
    numeric value. minimum density to consider CT.
}
  \item{minEvents}{
    %%     ~~Describe \code{minEvents} here~~
    numeric value. minimum events for the segmentation of breakpoints.
}
  \item{MC.CORES}{
    %%     ~~Describe \code{MC.CORES} here~~
    numeric value. number of cpu cores to use for computation.
}
}
\details{
  %%  ~~ If necessary, more details than the description above ~~
  Threshold values should depend on the platform/caller used. The default
  values are suited for WGS and Battenberg. 
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
This function returns a data.frame with all regions across samples
}
\references{
%% ~put references to the literature/web site here ~
The ICGC/TCGA Pan-Cancer Analysis of Whole Genomes Consortium. Pan-cancer analysis of whole genomes. Nature 578, 82–93 (2020). https://doi.org/10.1038/s41586-020-1969-6 
}
\author{
%%  ~~who you are~~
maxime.tarabichi@ulb.be; IRIBHM, Universite Libre de Bruxelles (ULB).
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~

}
\examples{

library(CTCallR)

## example run:
chromothripsis <-  CTCallR(CNA_file="./myCNA_calls.tsv",
                           SV_file="./mySV_calls.tsv",
                           outputfile="./myCT_calls.tsv")
head(chromothripsis)
head(colnames(chromothripsis))
table(chromothripsis$Calls)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
