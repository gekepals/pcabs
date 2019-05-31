\name{find_edge_between_groups}
\alias{find_edge_between_groups}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{find_edge_between_groups
}
\description{
Function to find the edges between two groups of variables.
}
\usage{
find_edge_between_groups(pdag, group1, group2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pdag}{
%%     ~~Describe \code{pdag} here~~
}
  \item{group1}{
%%     ~~Describe \code{group1} here~~
}
  \item{group2}{
%%     ~~Describe \code{group2} here~~
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
function (pdag, group1, group2)
{
    t = 0
    abs_group1 <- list()
    abs_group2 <- list()
    for (i in group1) {
        conn_list = find_conn_edges(pdag, i)
        for (z in conn_list) {
            if (z \%in\% group2) {
                t = t + 1
                abs_group1[t] <- i
                abs_group2[t] <- z
            }
        }
    }
    absList <- list(abs_group1 = abs_group1, abs_group2 = abs_group2)
    return(absList)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
