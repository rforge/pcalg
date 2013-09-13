\name{pdsep}
\alias{pdsep}
\title{Compute Possible-D-SEP for Each Node and Adapt Graph Accordingly}
\description{
  Computes the Possible-D-SEP for each node
  and tests conditional independence given all subsets of
  Possible-D-Sep. The graph is updated accordingly.
}
\usage{
pdsep(skel, suffStat, indepTest, p, sepset, alpha, pMax, m.max = Inf,
pdsep.max = Inf, NAdelete = TRUE, unfVect = NULL, biCC = FALSE, verbose = FALSE)
}
\arguments{
  \item{skel}{Graph object returned by \code{\link{skeleton}}.}
  \item{suffStat}{sufficient statistic: A list containing all necessary
    elements for making conditional independence decisions using
    function \code{indepTest}.}
  \item{indepTest}{predefined function for testing conditional independence. The
    function is internally called as \code{indepTest(x,y,S,suffStat)} for
    testing conditional independence of \code{x} and \code{y} given
    \code{S}. Here, \code{x} and \code{y} are node numbers of the
    adjacency matrix, \code{S} is a (possibly empty) vector of node
    numbers of the adjacency matrix and \code{suffStat} is a list
    containing all relevant elements for making conditional independence
    decisions. The return value of \code{indepTest} is the p-value of
    the test for conditional independence.}
  \item{p}{Number of variables.}
  \item{sepset}{List of length \code{p}; each element of the list
    contains another list of length \code{p}. The element
    \code{sepset[[x]][[y]]} contains the separation set that made the edge
    between \code{x} and \code{y} drop out. This object is thought to be
    obtained from a \code{pcAlgo}-object or \code{fciAlgo}-object.}
  \item{alpha}{Significance level for the individual conditional
    independence test.s}
  \item{pMax}{Matrix with the maximal p-values of conditional
    independence tests in a previous call of \code{\link{skeleton}},
    \code{\link{pc}} or \code{\link{fci}} which produced \code{G}. This
    object is thought to be obtained from a \code{pcAlgo}-object or
    \code{fciAlgo}-object.}
  \item{m.max}{Maximum size of the conditioning sets that are considered in the
    conditional independence tests.}
  \item{pdsep.max}{Maximum size of Possible-D-SEP for which subsets are
    considered as conditioning sets in the conditional independence
    tests. If the nodes \code{x} and \code{y} are adjacent in the graph
    and the size of Possible-D-SEP(\code{x})\ {\code{x},\code{y}},
    estimated normally or using biconnected components, see later, is
    bigger than pdsep.max, no conditional independence tests are
    performed between \code{x} and its neighbour \code{y} given subsets of
    Possible-D-SEP(\code{x}) and the edge is simply left in the graph.
    Note that if pdsep.max is less than Inf, the final PAG
    may be a supergraph than the one computed with pdsep.max = Inf,
    because less tests may have been performed in the former.}
  \item{NAdelete}{If indepTest returns \code{NA} and this option is
    \code{TRUE}, the corresponding edge is deleted. If this option is
    \code{FALSE}, the edge is not deleted.}
  \item{unfVect}{vector containing numbers that encode the unfaithful
    triple (as returned by \code{\link{pc.cons.intern}}). This is
    needed in the conservative FCI.}
  \item{biCC}{logical; if \code{TRUE}, only nodes on paths between nodes
    a and c are considered to be in sepset(a,c). Uses biconnected components,
    \code{\link[RBGL]{biConnComp}} from \pkg{RBGL}.}
  \item{verbose}{If \code{TRUE}, detailed output is provided.}
}
\value{A list with several elements:
  \item{G}{Updated adjacency matrix}
  \item{sepset}{Updated sepset}
  \item{pMax}{Updated matrix containing maximal p-values}
  \item{allPdsep}{Possible-D-Sep for each node}
  \item{max.ord}{Maximal order of conditioning sets during independence
    tests}
  \item{n.edgetests}{Number of conditional edgetests performed grouped
    by size of conditioning set.}
}

\details{For a given graph G, a node y is in Possible-D-SEP(x) if x!=y
  and there is an undirected path U between x and y in G such that for
  every subpath a,b,c of U either b is a collider on the subpath, or b
  is not a definite noncollider on U and a,b and c form a triangle in G.

  Each pair of nodes x and y which is connected by
  an edge is tested for conditional independence given every subset in
  Possible-D-SEP of x or of y. The conditional independence is tested on
  significance level alpha by using the test given in
  \code{indepTest}. If the pair of nodes is judged to be independent
  given set S, then S is recorded in the sepset of x,y and in the sepset
  of y,x and the edge is deleted. Otherwise, the edge is not deleted and
  no entry to the sepset is made.

  To make the code more efficient, tests that have been done during
  finding the skeleton are not done again.

  To make the estimation of Possible-D-SEP order-independent (because
  some edge deletion can influence sets computed later on), all sets are
  estimated at once and are not update after possible edge removals.}

\references{
  P. Spirtes, C. Glymour and R. Scheines (2000).
  \emph{Causation, Prediction, and Search}, 2nd edition. The MIT Press.
}
\seealso{\code{\link{qreach}} for a fast way of finding Possible-D-Sep
  for a given node; \code{\link{fci}} which uses \code{\link{pdsep}}.
}
\author{
  Markus Kalisch (\email{kalisch@stat.math.ethz.ch}) and Diego Colombo.
}
\examples{
p <- 10
## generate and draw random DAG:
set.seed(44)
myDAG <- randomDAG(p, prob = 0.2)

## generate 10000 samples of DAG using gaussian distribution
library(RBGL)
n <- 10000
d.mat <- rmvDAG(n, myDAG, errDist = "normal")

## estimate skeleton
indepTest <- gaussCItest
suffStat <- list(C = cor(d.mat), n = n)
alpha <- 0.01
skel <- skeleton(suffStat, indepTest, p, alpha)

## prepare input for pdsep
sepset <- skel@sepset
pMax <- skel@pMax

## call pdsep to find Possible-D-Sep and enhance the skeleton
pdsepRes <- pdsep(skel@graph, suffStat, indepTest, p, sepset, alpha,
pMax, verbose = TRUE)
## call pdsep with biconnected components to find Possible-D-Sep and enhance the skeleton
pdsepResBicc <- pdsep(skel@graph, suffStat, indepTest, p, sepset, alpha,
pMax, biCC= TRUE, verbose = TRUE)
}
