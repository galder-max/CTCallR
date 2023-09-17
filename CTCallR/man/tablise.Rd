\name{tablise}
\alias{tablise}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
tablise
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
tablise(cp.chr, i, nms)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cp.chr}{
%%     ~~Describe \code{cp.chr} here~~
}
  \item{i}{
%%     ~~Describe \code{i} here~~
}
  \item{nms}{
%%     ~~Describe \code{nms} here~~
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
function (cp.chr, i, nms) 
{
    .tablise <- function(cp.chr, i, nms, ww) {
        rbind(c(as.character(nms), as.character(cp.chr$start[1]), 
            as.character(cp.chr$end[1]), as.character(paste0("chr", 
                i, cp.chr$arms[1])), as.character(cp.chr$densities[1]), 
            as.character(cp.chr$Nevents[1]), as.character(cp.chr$testExpRegion[[1]]), 
            as.character(paste(cp.chr$cnStates[[ww[1]]]$states, 
                sep = "", collapse = ";")), as.character(paste(cp.chr$cnStates[[ww[1]]]$sizes, 
                sep = "", collapse = ";")), as.character(cp.chr$cnStates[[ww[1]]]$coveragemaxStates.modeSegs), 
            as.character(cp.chr$cnStates[[ww[1]]]$coveragemaxStates.maxCov)), 
            c(as.character(nms), as.character(cp.chr$start[2]), 
                as.character(cp.chr$end[2]), as.character(paste0("chr", 
                  i, cp.chr$arms[2])), as.character(cp.chr$densities[2]), 
                as.character(cp.chr$Nevents[2]), as.character(cp.chr$testExpRegion[[2]]), 
                as.character(paste(cp.chr$cnStates[[ww[2]]]$states, 
                  sep = "", collapse = ";")), as.character(paste(cp.chr$cnStates[[ww[2]]]$sizes, 
                  sep = "", collapse = ";")), as.character(cp.chr$cnStates[[ww[2]]]$coveragemaxStates.modeSegs), 
                as.character(cp.chr$cnStates[[ww[2]]]$coveragemaxStates.maxCov)))
    }
    ww <- which(!sapply(cp.chr$cnStates, is.null))[1:2]
    if (all(c("p", "q") \%in\% cp.chr$arms)) 
        return(.tablise(cp.chr, i, nms, ww))
    c(as.character(nms), as.character(cp.chr$start), as.character(cp.chr$end), 
        as.character(paste0("chr", i, cp.chr$arms)), as.character(cp.chr$densities), 
        as.character(cp.chr$Nevents), as.character(cp.chr$testExpRegion[[1]]), 
        as.character(paste(cp.chr$cnStates[[ww[1]]]$states, sep = "", 
            collapse = ";")), as.character(paste(cp.chr$cnStates[[ww[1]]]$sizes, 
            sep = "", collapse = ";")), as.character(cp.chr$cnStates[[ww[1]]]$coveragemaxStates.modeSegs), 
        as.character(cp.chr$cnStates[[ww[1]]]$coveragemaxStates.maxCov))
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
