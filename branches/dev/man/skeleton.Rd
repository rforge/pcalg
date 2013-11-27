\name{skeleton}
\alias{skeleton}
\title{Estimate the skeleton of a DAG without latent and selection
  variables using the PC Algorithm and estimate an initial skeleton of a
  DAG with arbitrarily many latent and selection variables using the FCI
  and the RFCI algorithms.}
\description{
  If used in the PC algorithm, it estimates the order-independent
  \dQuote{skeleton} of a directed acyclic graph (DAG) from observational
  data.  When used in the FCI and RFCI algorithms, this function estimates
  only an initial order-independent \dQuote{skeleton}.  Because of the
  presence of latent and selection variables, to find the final skeleton
  those algorithms need to perform additional tests later on and
  consequently some edges can be further deleted. 
}
\usage{
skeleton(suffStat, indepTest, p, alpha, fixedGaps = NULL,
         fixedEdges = NULL, NAdelete = TRUE, m.max = Inf, verbose = FALSE)
}
\arguments{
  \item{suffStat}{Sufficient statistics: List containing all necessary
    elements for the conditional independence decisions in the
    function \code{indepTest}.}
  \item{indepTest}{Predefined function for testing conditional independence. The
    function is internally called as \code{indepTest(x,y,S,suffStat)}, and
    tests conditional independence of \code{x} and \code{y} given
    \code{S}. Here, \code{x} and \code{y} are variables, and \code{S} is
    a (possibly empty) vector of variables (all variables are denoted
    by their column numbers
    in the adjacency matrix). \code{suffStat} is a list containing
    all relevant elements for the conditional independence
    decisions. The return value of \code{indepTest} is the p-value of
    the test for conditional independence.}
  \item{p}{Number of variables.}
  \item{alpha}{Significance level for the individual conditional
    independence tests.}
  \item{fixedGaps}{A logical matrix of dimension p*p. If entry
    \code{[i,j]} or \code{[j,i]} (or both) are TRUE, the edge i-j is
    removed before starting the algorithm. Therefore, this edge is
    guaranteed to be absent in the resulting graph.}
  \item{fixedEdges}{A logical matrix of dimension p*p. If entry
    \code{[i,j]} or \code{[j,i]} (or both) are TRUE, the edge i-j is
    never considered for removal. Therefore, this edge is
    guaranteed to be present in the resulting graph.}
  \item{NAdelete}{If indepTest returns \code{NA} and this option is
    \code{TRUE}, the corresponding edge is deleted. If this option is
    \code{FALSE}, the edge is not deleted.}
  \item{m.max}{Maximal size of the conditioning sets that are considered in the
    conditional independence tests.}
  \item{verbose}{If \code{TRUE}, detailed output is provided.}
}
\value{An object of \code{\link{class}} \code{"pcAlgo"} (see
  \code{\linkS4class{pcAlgo}}) containing an estimate of the skeleton of
  the underlying DAG, the conditioning sets that led to edge removals
  (sepset) and several other parameters.
}
\details{
  Under the assumption that the distribution of the observed variables
  is faithful to a DAG and that there are no latent and selection
  variables, this function estimates the skeleton of the DAG. The
  skeleton of a DAG is the undirected graph resulting from removing all
  arrowheads from the DAG. Edges in the skeleton of a DAG have the
  following interpretation: there is an edge between i and j if and only
  if variables i and j are conditionally dependent given S for all
  possible subsets S of the remaining nodes.

  Under the assumption that the distribution of the observed variables
  is faithful to a DAG with arbitrarily many latent and selection
  variables, this function estimates the initial skeleton of the
  DAG. Edges in this initial skeleton of a DAG have the
  following interpretation: there is an edge between i and j if and only
  if variables i and j are conditionally dependent given S for all
  possible subsets S of the neighbours of i and of the neighbours of j.

  The data are not required to follow a specific distribution,
  but one should make sure that the conditional indepedence test used in
  \code{indepTest} is appropriate for the data. Pre-programmed versions
  of \code{indepTest} are available for Gaussian data
  (\code{\link{gaussCItest}}), discrete data (\code{\link{disCItest}}),
  and binary data (see \code{\link{binCItest}}). Users can also specify
  their own \code{indepTest} function. 

  The PC algorithm (Spirtes, Glymour and Scheines, 2000)
  is known to be order-dependent, in the sense that the
  output depends on the order in which the variables are
  given. Therefore, Colombo and Maathuis (2013) proposed a simple
  modification, called PC-stable, that yields order-independent
  adjacencies in the skeleton. In this function we implement their
  modified algorithm (the old order-dependent implementation can be
  found in version 1.1-5).

  The algorithm starts with a complete undirected graph. In each
  step, it visits all pairs (\code{i},\code{j}) of adjacent nodes in the
  current graph, and determines based on conditional independence tests
  whether the edge i-j should be removed. In particular, for each step
  \code{m} (\code{m}=0,1,...) of the size of the conditioning sets, the
  algorithm at first determines the neighbours \code{a(i)} of each node
  \code{i} in the graph. Then, the algorithm visits all pairs (\code{i},
  \code{j}) of adjacent nodes in the current graph, and the edge between
  \code{i} and \code{j} is kept if and only if the null hypothesis
  "\code{i} and \code{j} are
  conditionally independent given S" is
  rejected at significance level \code{alpha} for all subsets \code{S} of size
  \code{m} of \code{a(i)} and of \code{a(j)} (as judged by the function
  \code{indepTest}). The difference with the original PC algorithm is
  that the neighborhoods \code{a(i)} are kept fixed within each value of
  \code{m}, and this makes the algorithm order-independent. The
  algorithm stops when \code{m} is larger than the largest neighbourhood
  size of all nodes, or when \code{m} has reached the limit \code{m.max}
  defined by the user. 

  Since the FCI (Spirtes, Glymour and Scheines, 2000) and RFCI (Colombo
  et al., 2012) algorithms are build up from the PC algorithm, they are also
  order-dependent in the skeleton. To resolve their order-dependence
  issues in the skeleton is more involved, see Colombo and Maathuis
  (2013). However, from the description above, this function estimates an
  initial order-independent skeleton in these algorithms (for additional
  details on how to make the final skeleton of FCI fully
  order-independent see \code{\link{fci}} and Colombo and Maathuis (2013)).   

  The information in \code{fixedGaps} and \code{fixedEdges} is used as follows.
  The gaps given in \code{fixedGaps} are introduced in the very beginning of
  the algorithm by removing the corresponding edges from the complete
  undirected graph. Pairs (\code{i},\code{j}) in \code{fixedEdges} are
  skipped in all steps of the algorithm, so that these edges remain in the
  graph.

  Note: Throughout, the algorithm works with the column positions of
  the variables in the adjacency matrix, and not with the names of
  the variables.
}
\references{
  D. Colombo and M.H. Maathuis (2013).
  \emph{Order-independent constraint-based causal structure learning},
  (arXiv:1211.3295v2)

  D. Colombo, M. H. Maathuis, M. Kalisch, T. S. Richardson
  (2012). Learning high-dimensional directed acyclic graphs with latent
  and selection variables. \emph{Ann. Statist.} \bold{40}, 294-321. 

  M. Kalisch and P. Buehlmann (2007).
  \emph{Estimating high-dimensional directed acyclic graphs with the
  PC-algorithm}, JMLR \bold{8} 613-636.
  
  P. Spirtes, C. Glymour and R. Scheines (2000).
  \emph{Causation, Prediction, and Search}, 2nd edition, MIT Press.     
}

\seealso{\code{\link{pc}} for generating a partially directed graph
  using the PC algorithm; \code{\link{fci}} for generating a partial
  ancestral graph using the FCI algorithm; \code{\link{rfci}} for
  generating a partial ancestral graph using the RFCI
  algorithm;\code{\link{gaussCItest}}, \code{\link{disCItest}},
  \code{\link{binCItest}} and \code{\link{dsepTest}} as examples for
  \code{indepTest}. 
}
\author{
  Markus Kalisch (\email{kalisch@stat.math.ethz.ch}), Martin Maechler,
  and Diego Colombo.
}
\examples{
##################################################
## Using Gaussian Data
##################################################
## Load predefined data
data(gmG)
n <- nrow(gmG$x)
p <- ncol(gmG$x)

## define independence test (partial correlations)
indepTest <- gaussCItest
## define sufficient statistics
suffStat <- list(C = cor(gmG$x), n = n)
## estimate Skeleton
alpha <- 0.01
skeleton.fit <- skeleton(suffStat, indepTest, p, alpha, verbose = TRUE)
if (require(Rgraphviz)) {
## show estimated Skeleton
par(mfrow=c(1,2))
plot(skeleton.fit, main = "Estimated Skeleton")
plot(gmG$g, main = "True DAG")
}
##################################################
## Using d-separation oracle
##################################################
## define independence test
indepTest <- dsepTest
## define sufficient statistics (d-separation oracle)
suffStat <- list(g = gmG$g, jp = RBGL::johnson.all.pairs.sp(gmG$g))
## estimate Skeleton
alpha <- 0.01 ## value is irrelevant as dsepTest returns either 0 or 1
fit <- skeleton(suffStat, indepTest, p, alpha)
if (require(Rgraphviz)) {
## show estimated Skeleton
plot(fit, main = "Estimated Skeleton")
plot(gmG$g, main = "True DAG")
}
##################################################
## Using discrete data
##################################################
## Load data
data(gmD)
p <- ncol(gmD$x)
## define independence test (G^2 statistics)
indepTest <- disCItest
## define sufficient statistics
suffStat <- list(dm = gmD$x, nlev = c(3,2,3,4,2), adaptDF = FALSE)
## estimate Skeleton
alpha <- 0.01
skeleton.fit <- skeleton(suffStat, indepTest, p, alpha, verbose = TRUE)
if (require(Rgraphviz)) {
## show estimated Skeleton
par(mfrow = c(1,2))
plot(skeleton.fit, main = "Estimated Skeleton")
plot(gmD$g, main = "True DAG")
}

##################################################
## Using binary data
##################################################
## Load binary data
data(gmB)
p <- ncol(gmB$x)
## define independence test
indepTest <- binCItest
## define sufficient statistics
suffStat <- list(dm = gmB$x, adaptDF = FALSE)
## estimate Skeleton
alpha <- 0.01
skeleton.fit <- skeleton(suffStat, indepTest, p, alpha, verbose = TRUE)
if (require(Rgraphviz)) {
## show estimated Skeleton
par(mfrow = c(1,2))
plot(skeleton.fit, main = "Estimated Skeleton")
plot(gmB$g, main = "True DAG")
}
}
\keyword{multivariate}
\keyword{models}
\keyword{graphs}
