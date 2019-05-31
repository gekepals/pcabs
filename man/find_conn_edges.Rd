\name{find_conn_edges}
\alias{find_conn_edges}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{find_conn_edges
}
\description{
Function to find the variables connected to a certain variable.
}
\usage{
find_conn_edges(pdag, row)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pdag}{
%%     ~~Describe \code{pdag} here~~
}
  \item{row}{
%%     ~~Describe \code{row} here~~
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
function (pdag, row)
{
    conn_list = list()
    t = 0
    z = 0
    for (i in pdag[row, ]) {
        t = t + 1
        if (i == 1) {
            z = z + 1
            conn_list[z] = colnames(pdag)[t]
        }
    }
    return(conn_list)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
