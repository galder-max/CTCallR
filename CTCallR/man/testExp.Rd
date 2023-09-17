\name{testExp}
\alias{testExp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
testExp
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
testExp(chrCN, N = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{chrCN}{
%%     ~~Describe \code{chrCN} here~~
}
  \item{N}{
%%     ~~Describe \code{N} here~~
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
function (chrCN, N = NULL) 
{
    require(vcd)
    require(MASS)
    sizes <- diff(chrCN[, "end"])
    if (length(sizes) < 4) 
        stop("number of breakpoints too low")
    if (any(sizes < 0)) {
        print(chrCN)
        stop("non-ordered chromosomal positions")
    }
    fit <- fitdistr(sizes, "exponential")
    ks.test(sizes, pexp, rate = fit$estimate)$p.value
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
