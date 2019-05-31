\name{find_pattern}
\alias{find_pattern}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{find_pattern
}
\description{
Function to produce the pdag from the skeleton.
}
\usage{
find_pattern(skel)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{skel}{
%%     ~~Describe \code{skel} here~~
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
function (skel)
{
    verbose <- TRUE
    res <- skel
    g <- as(skel, "matrix")
    p <- as.numeric(dim(g)[1])
    pdag <- g
    ind <- which(g == 1, arr.ind = TRUE)
    for (i in seq_len(nrow(ind))) {
        cat("\n", "i = ", i, ":")
        x <- ind[i, 1]
        cat("\n", x)
        y <- ind[i, 2]
        cat("\n", y, "\n")
        allz <- setdiff(which(g[y, ] == 1), x)
        cat("z:", allz, "\n")
        for (z in allz) {
            if (g[x, z] == 0 && !(y \%in\% skel@sepset[[x]][[z]] ||
                y \%in\% skel@sepset[[z]][[x]])) {
                if (verbose) {
                  cat("\n", x, "->", y, "<-", z, "\n")
                  cat("Sxz=", skel@sepset[[z]][[x]], "Szx=",
                    skel@sepset[[x]][[z]])
                }
                pdag[x, y] <- pdag[z, y] <- 1
                pdag[y, x] <- pdag[y, z] <- 0
            }
        }
    }
    return(pdag)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
