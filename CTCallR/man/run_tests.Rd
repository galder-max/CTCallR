\name{run_tests}
\alias{run_tests}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
run_tests
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
run_tests(ct, lSV, samples)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ct}{
%%     ~~Describe \code{ct} here~~
}
  \item{lSV}{
%%     ~~Describe \code{lSV} here~~
}
  \item{samples}{
%%     ~~Describe \code{samples} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (ct, lSV, samples) 
{
    chr <- gsub("chr", "", gsub("p", "", gsub("q", "", ct[4])))
    if (chr == "23") 
        chr <- "X"
    svs <- lSV[[which(samples == ct["samplename"])]]
    require(GenomicRanges)
    grSV1 <- GRanges(svs[, 2], IRanges(svs[, 4], svs[, 4]))
    grSV2 <- GRanges(svs[, 5], IRanges(svs[, 7], svs[, 7]))
    grCT <- GRanges(chr, IRanges(as.numeric(ct[2]), as.numeric(ct[3])))
    ov1 <- findOverlaps(grSV1, grCT)
    ov2 <- findOverlaps(grSV1, grCT)
    keep <- unique(c(queryHits(ov1), queryHits(ov2)))
    svs <- svs[keep, ]
    list(testclass(svs), testorder(svs[svs[, 2] == chr & svs[, 
        5] == chr, ]), length(keep))
  }
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
