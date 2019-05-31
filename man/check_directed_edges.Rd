\name{check_directed_edges}
\alias{check_directed_edges}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{check_directed_edges
}
\description{
Function to check for directed edges between two groups of variables.
}
\usage{
check_directed_edges(pdag, abs_group1, abs_group2)
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
    ind <- which((pdag == 1 & t(pdag) == 0), arr.ind = TRUE)
    dir_edge <- FALSE
    ind_no <- list()
    t <- 1
    for (i in seq_len(nrow(ind))) {
        if (ind[i, 1] \%in\% abs_group1) {
            if (ind[i, 2] \%in\% abs_group2) {
                cat("directed edge to the right found on t = ",
                  i, "\n")
                dir_edge = TRUE
                ind_no[[t]] <- ind[i, 1]
                t <- t + 1
                ind_no[[t]] <- ind[i, 2]
                t <- t + 1
            }
        }
        else if (ind[i, 2] \%in\% abs_group1) {
            if (ind[i, 1] \%in\% abs_group2) {
                cat("directed edge to the left found on t = ",
                  i, "\n")
                dir_edge = TRUE
                ind_no[[t]] <- ind[i, 1]
                t <- t + 1
                ind_no[[t]] <- ind[i, 2]
                t <- t + 1
            }
        }
    }
    if (dir_edge) {
        return(ind_no)
    }
    else {
        return(NULL)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
