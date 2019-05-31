\name{is_direction_possible}
\alias{is_direction_possible}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{is_direction_possible
}
\description{
This function checks whether a certain direction of the edges is possible.
It first check whether the direction alone causes any more v-structures.
Then, it applies the MEC-Rules of the pcalg algorithm, and checks whether
this causes any more v-structures.
If both don't cause more v-structures, the direction is allowed.
}
\usage{
is_direction_possible(pdag, v_structures)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pdag}{
%%     ~~Describe \code{pdag} here~~
}
  \item{v_structures}{
%%     ~~Describe \code{v_structures} here~~
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
function (pdag, v_structures)
{
    possible_direction <- TRUE
    pdag_v <- number_v_nodes(pdag)
    if (v_structures != pdag_v) {
        possible_direction <- FALSE
    }
    if (possible_direction) {
        pdag <- apply_mec_rules(pdag)
        pdag_v <- number_v_nodes(pdag)
        if (v_structures != pdag_v) {
            possible_direction <- FALSE
        }
    }
    return(possible_direction)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
