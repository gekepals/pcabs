\name{add_abstraction}
\alias{add_abstraction}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{add_abstraction
}
\description{
Directs the edges between a pair of clusters.
There are 3 options when an abstraction is added:
1) there is no directed edge in the abstraction. This means both directions must be tried.
2) there is 1 directed edge in the abstraction. This means all edges must be pointed that way.
3) there are more than 1 directed edges in the abstraction. This means that first, it must
   be checked whether they all point in the same direction. Then, it must be checked if there
   are undirected edges left. If so, they must point in the same direction as well.
}
\usage{
add_abstraction(pdag, abs_group1, abs_group2)
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
    abs_group1 <- as.numeric(noquote(gsub("[^0-9]", "", abs_group1)))
    abs_group2 <- as.numeric(noquote(gsub("[^0-9]", "", abs_group2)))
    dir_list <- check_directed_edges(pdag, abs_group1, abs_group2)
    v_structures <- number_v_nodes(pdag)
    if (length(dir_list) == 0) {
        print("no directed edges found on the abstraction")
        pdag1 <- change_direction(pdag, abs_group1, abs_group2)
        dir1 <- is_direction_possible(pdag1, v_structures)
        pdag2 <- change_direction(pdag, abs_group2, abs_group1)
        dir2 <- is_direction_possible(pdag2, v_structures)
        if (dir1 && dir2) {
            cat("Two possible directions of the edges. Need to check for any necessary edges among them.")
            end_pdag_list <- list()
            col <- ncol(pdag)
            row <- nrow(pdag)
            for (i in 1:length(pdag)) {
                if (pdag1[i] == 0 && pdag2[i] == 0) {
                  end_pdag_list[i] = 0
                }
                else {
                  end_pdag_list[i] = 1
                }
            }
            end_pdag <- matrix(end_pdag_list, nrow = row, ncol = col)
            return(end_pdag)
        }
        if (dir1 && !dir2) {
            cat("The direction from abs_group1 -> abs_group2 is the only possible direction.")
            pdag <- apply_mec_rules(pdag1)
            return(pdag)
        }
        if (!dir1 && dir2) {
            cat("The direction from abs_group1 <- abs_group2 is the only possible direction.")
            pdag <- apply_mec_rules(pdag2)
            return(pdag)
        }
        if (!dir1 && !dir2) {
            cat("Both the direction are not possible...")
            return(pdag)
        }
    }
    else if (length(dir_list) == 2) {
        print("one directed edge found. Need to convert other edges the same way.")
        group1_check <- check_direction(pdag, dir_list, abs_group1)
        if (!(group1_check)) {
            q <- abs_group1
            abs_group1 <- abs_group2
            abs_group2 <- q
        }
        pdag <- change_direction(pdag, abs_group1, abs_group2)
        pdag <- apply_mec_rules(pdag)
        return(pdag)
    }
    else {
        print("multiple directed edges found. Need to check if they match.")
        leftrightlist <- list_dir_edges(dir_list)
        result <- match_dir_edges(leftrightlist$left, leftrightlist$right,
            abs_group1)
        if (result == "error") {
            cat("Error found!")
        }
        else if (result == "left") {
            cat("the directed edges match with the pattern! Leftlist is in absgroup1.")
            for (i in 1:abs_group1) {
                if (!(i \%in\% leftrightlist$left)) {
                  cat("Not all edges are directed yet.")
                  group1_check <- check_direction(pdag, dir_list,
                    abs_group1)
                  if (!(group1_check)) {
                    q <- abs_group1
                    abs_group1 <- abs_group2
                    abs_group2 <- q
                  }
                  pdag <- change_direction(pdag, abs_group1,
                    abs_group2)
                  return(pdag)
                }
                else {
                  cat("all edges are already directed equally. The abstraction doesn't give extra information.")
                  return(pdag)
                }
            }
        }
        else if (result == "right") {
            cat("the directed edges match with the pattern! Rightlist is in absgroup1.")
            for (i in 1:abs_group1) {
                if (!(i \%in\% leftrightlist$right)) {
                  cat("Not all edges are directed yet.")
                  group1_check <- check_direction(pdag, dir_list,
                    abs_group1)
                  if (!(group1_check)) {
                    q <- abs_group1
                    abs_group1 <- abs_group2
                    abs_group2 <- q
                  }
                  pdag <- change_direction(pdag, abs_group1,
                    abs_group2)
                  return(pdag)
                }
                else {
                  cat("all edges are already directed equally. The abstraction doesn't give extra information.")
                  return(pdag)
                }
            }
        }
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
