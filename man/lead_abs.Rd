\name{lead_abs}
\alias{lead_abs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{lead_abs
}
\description{
Takes as input the pdag and all the clusters (all Zs).
Loops over the groups to decide which are pairs (i.e. which have connected edges).
Decides the order of directing edges between pairs.
First, the pairs are directed that have 1 or more directed edges.
Last, all pairs that have no directed edges are directed.
All pairs are sent to the add_abstraction function to direct the edges.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
lead_abs(pdag, abs_groups)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pdag}{
%%     ~~Describe \code{pdag} here~~
}
  \item{abs_groups}{
%%     ~~Describe \code{abs_groups} here~~
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
function (pdag, abs_groups)
{
    undirected_list <- list()
    for (i in 1:length(abs_groups)) {
        cat("i: ", i)
        for (z in i:length(abs_groups)) {
            if ((z + 1) <= length(abs_groups)) {
                absList <- find_edge_between_groups(pdag, abs_groups[[i]],
                  abs_groups[[z + 1]])
                if (length(absList$abs_group1) > 0) {
                  group1 <- noquote(gsub("[^0-9]", "", absList$abs_group1))
                  group2 <- noquote(gsub("[^0-9]", "", absList$abs_group2))
                  dir_list <- check_directed_edges(pdag, group1,
                    group2)
                  if (!is.null(dir_list)) {
                    print("directed edges found between groups")
                    pdag <- add_abstraction(pdag, absList$abs_group1,
                      absList$abs_group2)
                  }
                  else {
                    print("no directed edges found between groups")
                    undirected_list <- append(undirected_list,
                      absList)
                  }
                }
            }
        }
    }
    if (length(undirected_list) > 0) {
        print("start with undirected list")
        for (i in seq(1, length(undirected_list), by = 2)) {
            pdag <- add_abstraction(pdag, undirected_list[i]$abs_group1,
                undirected_list[i + 1]$abs_group2)
        }
    }
    return(pdag)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
