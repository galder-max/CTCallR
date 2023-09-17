\name{treatPCF}
\alias{treatPCF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
treatPCF
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
treatPCF(pcf, flagDensity = 1.5e+08/50, minEvents = 5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pcf}{
%%     ~~Describe \code{pcf} here~~
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
function (pcf, flagDensity = 1.5e+08/50, minEvents = 5) 
{
    npcf <- pcf
    whichFlaggedP <- which(pcf$mean < flagDensity & pcf$n.probes >= 
        minEvents & pcf$arm == "p")
    if (length(whichFlaggedP) > 0) {
        whichFlaggedP <- min(whichFlaggedP):max(whichFlaggedP)
        if (length(whichFlaggedP) > 1) 
            npcf <- npcf[-c(whichFlaggedP[2:length(whichFlaggedP)]), 
                ]
        i <- whichFlaggedP[1]
        npcf[i, "start.pos"] <- pcf[min(whichFlaggedP), "start.pos"]
        npcf[i, "end.pos"] <- pcf[max(whichFlaggedP), "end.pos"]
        npcf[i, "n.probes"] <- sum(pcf[whichFlaggedP, "n.probes"])
        wfp <- whichFlaggedP
        npcf[i, "mean"] <- sum(pcf[wfp, "mean"]/1000 * pcf[wfp, 
            "n.probes"])/sum(pcf[wfp, "n.probes"]) * 1000
    }
    pcf <- npcf
    whichFlaggedP <- which(pcf$mean < flagDensity & pcf$n.probes >= 
        minEvents & pcf$arm == "q")
    if (length(whichFlaggedP) > 0) {
        whichFlaggedP <- min(whichFlaggedP):max(whichFlaggedP)
        if (length(whichFlaggedP) > 1) 
            npcf <- npcf[-c(whichFlaggedP[2:length(whichFlaggedP)]), 
                ]
        i <- whichFlaggedP[1]
        npcf[i, "start.pos"] <- pcf[min(whichFlaggedP), "start.pos"]
        npcf[i, "end.pos"] <- pcf[max(whichFlaggedP), "end.pos"]
        npcf[i, "n.probes"] <- sum(pcf[whichFlaggedP, "n.probes"])
        wfp <- whichFlaggedP
        npcf[i, "mean"] <- sum(pcf[wfp, "mean"]/1000 * pcf[wfp, 
            "n.probes"])/sum(pcf[wfp, "n.probes"]) * 1000
    }
    return(npcf)
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
