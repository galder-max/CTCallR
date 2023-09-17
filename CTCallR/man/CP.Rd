\name{CP}
\alias{CP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
CP
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
CP(chrCN, flagDensity, minEvents)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{chrCN}{
%%     ~~Describe \code{chrCN} here~~
}
  \item{flagDensity}{
%%     ~~Describe \code{flagDensity} here~~
}
  \item{minEvents}{
%%     ~~Describe \code{minEvents} here~~
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
function (chrCN, flagDensity, minEvents) 
{
    if (nrow(chrCN) > minEvents) {
        pcf <- doPCF(chrCN, minEvents = minEvents)
        pcf <- treatPCF(pcf, flagDensity, minEvents)
        starts <- cumsum(c(1, pcf$n.probes[-c(nrow(pcf))]))
        ends <- cumsum(pcf$n.probes)
        whichFlagged <- which(pcf$mean < flagDensity & pcf$n.probes >= 
            minEvents)
        tE <- mytry(testExp(chrCN, N = NULL))
        tERegion <- lapply(whichFlagged, function(x) {
            res <- try(testExp(chrCN[starts[x]:ends[x], ], N = N), 
                silent = T)
            if (inherits(res, "try-error")) 
                return(NULL)
            res
        })
        flaggedtE <- tE < 0.05
        flaggedDensity <- length(whichFlagged) > 0
        flagged <- flaggedtE | flaggedDensity
        copyNstates <- getCNstates(pcf, whichFlagged, chrCN)
        return(list(starts = chrCN[starts[whichFlagged], "start"], 
            ends = chrCN[ends[whichFlagged], "end"], arms = pcf[whichFlagged, 
                "arm"], densities = pcf$mean[whichFlagged], Nevents = pcf$n.probes[whichFlagged], 
            cnStates = copyNstates, testExpChromosome = tE, testExpRegion = tERegion, 
            flaggedtE = flaggedtE, flaggedDensity = flaggedDensity, 
            flagged = flagged))
    }
    else return(list(flagged = FALSE, minEventsReached = FALSE))
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
