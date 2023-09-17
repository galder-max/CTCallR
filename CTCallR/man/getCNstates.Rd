\name{getCNstates}
\alias{getCNstates}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
getCNstates
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getCNstates(pcf, whichFlagged, chrCN, maxStates = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pcf}{
%%     ~~Describe \code{pcf} here~~
}
  \item{whichFlagged}{
%%     ~~Describe \code{whichFlagged} here~~
}
  \item{chrCN}{
%%     ~~Describe \code{chrCN} here~~
}
  \item{maxStates}{
%%     ~~Describe \code{maxStates} here~~
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
function (pcf, whichFlagged, chrCN, maxStates = 3) 
{
    starts <- c(0, cumsum(pcf$n.probes[-c(nrow(pcf))])) + 1
    ends <- c(cumsum(pcf$n.probes))
    nStates <- list()
    for (i in whichFlagged) {
        wF <- as.logical(1:nrow(chrCN) * 0)
        wF[starts[i]:ends[i]] <- TRUE
        states <- paste(chrCN[wF, "nMaj1"], chrCN[wF, "nMin1"], 
            sep = "-")
        sizes <- chrCN[wF, "end"] - chrCN[wF, "start"]
        isSubclonal <- chrCN[wF, "frac1"]
        ssizes <- tapply(sizes, as.factor(states), sum)/1000
        nsegs <- tapply(sizes, as.factor(states), length)
        fsizes <- ssizes/sum(ssizes)
        nS <- length(ssizes)
        if (nS < maxStates) 
            maxStates <- nS
        nStates[[i]] <- list(states = states, sizes = sizes, 
            ssizes = ssizes, isSubclonal = isSubclonal, subclonalFrac = sum(sizes/100 * 
                as.numeric(isSubclonal != 1), na.rm = T)/sum(sizes/100, 
                na.rm = T), fsizes = fsizes, nsegs = nsegs, coveragemaxStates.modeSegs = sum(fsizes[order(nsegs, 
                decreasing = T)][1:maxStates]), coveragemaxStates.maxCov = sum(sort(fsizes, 
                decreasing = T)[1:maxStates]))
    }
    nStates
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
