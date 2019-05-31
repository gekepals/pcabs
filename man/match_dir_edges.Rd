\name{match_dir_edges}
\alias{match_dir_edges}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{match_dir_edges
}
\description{
Function to check if variables in leftlist are in either absgroup1 or absgroup2.
This is necessary, otherwise there has been made a mistake in the given absgroups.
}
\usage{
match_dir_edges(leftlist, rightlist, absgroup1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{leftlist}{
%%     ~~Describe \code{leftlist} here~~
}
  \item{rightlist}{
%%     ~~Describe \code{rightlist} here~~
}
  \item{absgroup1}{
%%     ~~Describe \code{absgroup1} here~~
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
function (leftlist, rightlist, absgroup1)
{
    result <- "error"
    group1 <- FALSE
    group2 <- FALSE
    for (i in leftlist) {
        if (i \%in\% absgroup1) {
            group1 <- TRUE
        }
        else {
            group2 <- TRUE
        }
    }
    if (group1 && group2) {
        return(result)
    }
    if (group1) {
        for (i in rightlist) {
            if (i \%in\% absgroup1) {
                group1 <- FALSE
            }
            else {
                group2 <- TRUE
            }
        }
        if (group1 && group2) {
            cat("Nice! The abstraction matches with the pattern.")
            result <- "left"
            return(result)
        }
    }
    else if (group2) {
        for (i in rightlist) {
            if (i \%in\% absgroup1) {
                group1 <- TRUE
            }
            else {
                group2 <- FALSE
            }
        }
        if (group1 && group2) {
            cat("Nice! The abstraction matches witht he pattern.")
            result <- "right"
            return(result)
        }
        else {
            return(result)
        }
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
