\name{list_dir_edges}
\alias{list_dir_edges}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{list_dir_edges
}
\description{
Function to put the results of check_directed_edges function into two lists:
One list with all the nodes found on the left side, one with nodes on the right side.
The functions thus splits up the list found in check_directed_edges.
}
\usage{
list_dir_edges(dir_list)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dir_list}{
%%     ~~Describe \code{dir_list} here~~
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
function (dir_list)
{
    left_list <- list()
    right_list <- list()
    l <- 1
    r <- 1
    for (i in seq(1, length(dir_list), 2)) {
        left_list[l] <- dir_list[[i]]
        l <- l + 1
    }
    for (i in seq(2, length(dir_list), 2)) {
        right_list[r] <- dir_list[[i]]
        r <- r + 1
    }
    leftright_list <- list(left = left_list, right = right_list)
    return(leftright_list)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
