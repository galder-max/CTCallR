\name{testclass}
\alias{testclass}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
testclass
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
testclass(svs)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{svs}{
%%     ~~Describe \code{svs} here~~
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
function (svs) 
{
    require(EMT)
    classes <- svs[, "svclass"]
    types <- c("Deletion", "Head-to-head_inversion", "Tail-to-Tail_inversion", 
        "Tandem_duplication")
    observed <- as.vector(table(classes)[types])
    observed[is.na(observed)] <- 0
    Chi2 <- try(chisq.test(observed, p = rep(1/4, 4))$p.value, 
        silent = T)
    if (inherits(Chi2, "try-error")) 
        Chi2 <- NA
    Chi2
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
