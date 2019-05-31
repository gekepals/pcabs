\name{apply_mec_rules}
\alias{apply_mec_rules}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{apply_mec_rules
}
\description{
Function that applies the 3 rules of PC algorithm.
}
\usage{
apply_mec_rules(pdag)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pdag}{
%%     ~~Describe \code{pdag} here~~
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
function (pdag)
{
    verbose <- TRUE
    g <- as(skel, "matrix")
    p <- as.numeric(dim(g)[1])
    old_pdag <- matrix(0, p, p)
    while (!all(old_pdag == pdag)) {
        old_pdag <- pdag
        ind <- which((pdag == 1 & t(pdag) == 0), arr.ind = TRUE)
        cat("ind in Rule 1: ", ind)
        for (i in seq_len(nrow(ind))) {
            cat("i ", i)
            a <- ind[i, 1]
            cat("a ", a)
            b <- ind[i, 2]
            cat("b ", b)
            indC <- which((pdag[b, ] == 1 & pdag[, b] == 1) &
                (pdag[a, ] == 0 & pdag[, a] == 0))
            cat("indC in Rule 1: ", indC)
            if (length(indC) > 0) {
                pdag[b, indC] <- 1
                pdag[indC, b] <- 0
                if (verbose)
                  cat("\nRule 1:", a, "->", b, " and", b, "-",
                    indC, " where ", a, " and ", indC, " not connected : ",
                    b, "->", indC, "\n")
            }
        }
        ind <- which((pdag == 1 & t(pdag) == 1), arr.ind = TRUE)
        for (i in seq_len(nrow(ind))) {
            a <- ind[i, 1]
            b <- ind[i, 2]
            indC <- which((pdag[a, ] == 1 & pdag[, a] == 0) &
                (pdag[, b] == 1 & pdag[b, ] == 0))
            if (length(indC) > 0) {
                pdag[a, b] <- 1
                pdag[b, a] <- 0
                if (verbose)
                  cat("\nRule 2: Kette ", a, "->", indC, "->",
                    b, ":", a, "->", b, "\n")
            }
        }
        ind <- which((pdag == 1 & t(pdag) == 1), arr.ind = TRUE)
        for (i in seq_len(nrow(ind))) {
            a <- ind[i, 1]
            b <- ind[i, 2]
            indC <- which((pdag[a, ] == 1 & pdag[, a] == 1) &
                (pdag[, b] == 1 & pdag[b, ] == 0))
            if (length(indC) >= 2) {
                g2 <- pdag[indC, indC]
                if (length(g2) <= 1) {
                  g2 <- 0
                }
                else {
                  diag(g2) <- rep(1, length(indC))
                }
                if (any(g2 == 0)) {
                  pdag[a, b] <- 1
                  pdag[b, a] <- 0
                  if (verbose)
                    cat("\nRule 3:", a, "->", b, "\n")
                }
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
