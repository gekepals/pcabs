\name{change_direction}
\alias{change_direction}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{change_direction
}
\description{
Function that changes all undirected edges into directed edges.
The edges are directed as abs_group1 --> abs_group2.
}
\usage{
change_direction(pdag, abs_group1, abs_group2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pdag}{
%%     ~~Describe \code{pdag} here~~
}
  \item{abs_group1}{
%%     ~~Describe \code{abs_group1} here~~
}
  \item{abs_group2}{
%%     ~~Describe \code{abs_group2} here~~
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
function (pdag, abs_group1, abs_group2)
{
    for (i in seq_along(abs_group1)) {
        a <- abs_group1[i]
        b <- abs_group2[i]
        if (pdag[a, b] == 1 & pdag[b, a] == 1) {
            pdag[b, a] <- 0
        }
    }
    return(pdag)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
