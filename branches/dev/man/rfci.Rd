\name{rfci}
\alias{rfci}
\title{Estimate the equivalence class of a MAG (PAG) using the RFCI Algorithm}
\description{
  Estimate the equivalence class of a maximal ancestral graph (MAG) from
  observational data, using the RFCI-algorithm.
}

\usage{
rfci(suffStat, indepTest, p, alpha, verbose = FALSE, fixedGaps = NULL,
    fixedEdges = NULL, NAdelete = TRUE, m.max = Inf, rules = rep(TRUE,
10), conservative = FALSE, maj.rule = FALSE, cons.rules = FALSE, labels = NA)
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
  \item{verbose}{If \code{TRUE}, detailed output is provided.}
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
  \item{m.max}{Maximum size of the conditioning sets that are considered in the
    conditional independence tests.}
  \item{rules}{Logical vector of length 10 indicating which rules
    should be used when directing edges. The order of the rules is taken
    from Zhang (2009).}
  \item{conservative}{If \code{TRUE}, the
  triples are checked for ambiguity after the skeleton has been found.}
  \item{maj.rule}{If \code{TRUE}, the
  triples are checked for ambiguity after the skeleton has been found
  using the majority rule idea, which
  is less strict than the standard conservative. Note that to use the
  majority rule idea, the user has to set conservative = TRUE as well.} 
  \item{cons.rules}{If \code{TRUE}, an orientation rule that needs information
    on definite non-colliders is only applied, if the corresponding
    subgraph relevant for the rule does not involve an ambiguous triple.}
  \item{labels}{Vector of length p with node names. If \code{NA},
    \code{as.character(1:p)} is used instead.}
}
\value{An object of \code{\link{class}} \code{fciAlgo} (see
 \code{\linkS4class{fciAlgo}}) containing the estimated graph
 (in the form of an adjacency matrix with various possible edge marks),
 the conditioning sets that lead to edge removals (sepset) and several other
 parameters.
}

\details{
  This function is rather similar to \link{fci}. However, it does not
  compute any Possible-D-SEP sets and thus does not make tests
  conditioning on them. This makes RFCI much faster than FCI. The
  orientation rules for v-structures and rule 4 were modified in order
  to produce a PAG which, in the oracle version, is guaranteed to have
  the correct ancestral relationships.

  The argument \code{conservative} is used to estimate the conservative
  RFCI. If \code{conservative = TRUE} it is computed calling
  \code{pc.cons.intern} with option \code{version.unf = c(1,1)} after
  the initial skeleton has been computed and the algorithm is performing
  the additional local tests on the unshielded triples. The resulting
  information on sepset and unambiguous triples will be used for the
  subsequent orientation rules. Instead of using the conservative
  version, which is quite strict towards the v-structures,
  \cite{pc.indep} introduced a less strict version for the v-structures
  called majority rule, see \code{pc.cons.intern} for more details. This
  adaptation can be called using \code{conservative = TRUE} and
  \code{maj.rule = TRUE}. The purpose of this conservative versions is
  used in the hope to obtain better edge orientations.

  In contrast to the FCI algorithm,
  none of the order-dependence issues are fully resolved here, see
  Colombo and Maathuis (2013).}

\references{
  D. Colombo and M.H. Maathuis (2013). Order-independent
  constraint-based causal structure learning. arXiv preprint
  arXiv:1211.3295v1.
  
  D. Colombo, M. H. Maathuis, M. Kalisch, T. S. Richardson
  (2012). Learning high-dimensional directed acyclic graphs with latent
  and selection variables. \emph{Ann. Statist.} \bold{40}, 294-321. 
}
\seealso{\code{\link{fci}} for estimating a PAG using the FCI algorithm;
  \code{\link{skeleton}} for estimating a skeleton
  using the PC algorithm; \code{\link{pc}} for estimating a CPDAG using
  the PC algorithm; \code{\link{gaussCItest}},
  \code{\link{disCItest}}, \code{\link{binCItest}} and
  \code{\link{dsepTest}} as examples for \code{indepTest}.
}
\author{
  Diego Colombo and Markus Kalisch (\email{kalisch@stat.math.ethz.ch}).
}
\examples{
##################################################
## Example without latent variables
##################################################
set.seed(42)
p <- 7
## generate and draw random DAG :
myDAG <- randomDAG(p, prob = 0.4)

## find skeleton and PAG using the RFCI algorithm
suffStat <- list(C = cov2cor(trueCov(myDAG)), n = 10^9)
indepTest <- gaussCItest
res <- rfci(suffStat, indepTest, p, alpha = 0.9999, verbose=TRUE)


##################################################
## Example with hidden variables
## Zhang (2008), Fig. 6, p.1882
##################################################

## create the graph
p <- 4
amat <- t(matrix(c(0,1,0,0,1, 0,0,1,0,0, 0,0,0,1,0, 0,0,0,0,0, 0,0,0,1,0),5,5))
colnames(amat) <- rownames(amat) <- as.character(1:5)
L <- 1
V <- as.character(1:5)
edL <- vector("list",length=5)
names(edL) <- V
edL[[1]] <- list(edges=c(2,4),weights=c(1,1))
edL[[2]] <- list(edges=3,weights=c(1))
edL[[3]] <- list(edges=5,weights=c(1))
edL[[4]] <- list(edges=5,weights=c(1))
g <- new("graphNEL", nodes=V, edgeL=edL, edgemode="directed")

## compute the true covariance matrix of g
cov.mat <- trueCov(g)

## delete rows and columns belonging to latent variable L
true.cov <- cov.mat[-L,-L]

## transform covariance matrix into a correlation matrix
true.corr <- cov2cor(true.cov)

## find PAG with RFCI algorithm
## as dependence "oracle", we use the true correlation matrix in the
## function gaussCItest with a large "virtual sample size" and a large
## alpha
suffStat <- list(C = true.corr, n = 10^9)
indepTest <- gaussCItest
rfci.pag <- rfci(suffStat, indepTest, p, alpha = 0.9999, verbose=TRUE)

## define PAG given in Zhang (2008), Fig. 6, p.1882
corr.pag <- matrix(0,4,4)
corr.pag[1,] <- c(0,1,1,0)
corr.pag[2,] <- c(1,0,0,2)
corr.pag[3,] <- c(1,0,0,2)
corr.pag[4,] <- c(0,3,3,0)

## check if estimated and correct PAG are in agreement
all(corr.pag == rfci.pag@amat)

}
\keyword{multivariate}
\keyword{models}
\keyword{graphs}
