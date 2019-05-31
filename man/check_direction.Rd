\name{check_direction}
\alias{check_direction}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{check_direction
}
\description{
Function to check the direction of the variables in dir_list.
Assumptions: the variables in the abs_groups are directed equally,
and there is at least 1 directed edge in abs_group1.
(this has been checked before the function is called)
}
\usage{
check_direction(pdag, dir_list, abs_group1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pdag}{
%%     ~~Describe \code{pdag} here~~
}
  \item{dir_list}{
%%     ~~Describe \code{dir_list} here~~
}
  \item{abs_group1}{
%%     ~~Describe \code{abs_group1} here~~
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
function (pdag, dir_list, abs_group1)
{
    a <- dir_list[[1]]
    b <- dir_list[[2]]
    if (pdag[a, b] == 1 & pdag[b, a] == 0) {
        cat("the edge is directed as ", a, "->", b)
        if (a \%in\% abs_group1) {
            return(group1 <- TRUE)
        }
        else {
            return(group1 <- FALSE)
        }
    }
    if (pdag[b, a] == 1 & pdag[a, b] == 0) {
        cat("the edge is directed as ", a, "<-", b)
    }
    if (b \%in\% abs_group1) {
        return(group1 <- FALSE)
    }
    else {
        return(group1 <- TRUE)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
