\name{udag2apag}
\alias{udag2apag}
\title{Extend a partially oriented skeleton into an RFCI-PAG}
\description{
This function extends a skeleton in which the v-structures have been
already oriented into an RFCI Partial Ancestral Graph (PAG). While
orienting the edges, this function performs some additional conditional
independence tests in rule 4 to ensure correctness of the ancestral
relationships, and consequently some edges can be deleted. The result is
the final adjacency matrix indicating also the edge marks and the
updated sepsets.  
}
\usage{
udag2apag(apag, suffStat, indepTest, alpha, sepset, rules = rep(TRUE,
10), unfVect = NULL, verbose = FALSE)
}
\arguments{
  \item{apag}{adjacency matrix of the pcAlgo-object of size p*p where
    the v-structures have already been oriented.}
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
  \item{alpha}{Significance level for the individual conditional
    independence tests.}
  \item{sepset}{List of length p; each element of the list
    contains another list of length p. The element
    \code{sepset[[x]][[y]]} contains the separation set that made the edge
    between \code{x} and \code{y} drop out. This object is thought to be
    obtained from a \code{pcAlgo}-object.}
  \item{rules}{array of length 10 containing \code{TRUE} or \code{FALSE} for each
    rule. \code{TRUE} in position i means that rule i (Ri) will be used. Per
    default, all rules are set to \code{TRUE}.}
  \item{unfVect}{Vector containing numbers that encode the unfaithful
    triple. This is needed in the conservative and in the majority rule RFCI.}
  \item{verbose}{If \code{TRUE}, detailed output is provided.}
}
\details{
  The partially oriented skeleton (v-structures have been already
  oriented) is extended to an RFCI-PAG using adapted rules by Zhang (see
  References). This function is similar to \link{udag2pag} used to
  orient the skeleton into a PAG in the FCI algorithm. However, it is
  slightly more complicated because while orienting rule 4 additional
  conditional independence tests are performed to ensure correctness of
  the ancestral relationships and some edges can be deleted, see Colombo
  et al. (2012). To this purpose, \code{suffStat}, \code{indepTest}, and
  \code{alpha} are used as inputs. Note that because some edges can be
  deleted, the input adjacency matrix \code{apag} and the input
  separating sets \code{sepset} can change in the run of the algorithm.

  If \code{unfVect = NULL} the ten orientation rules are applied to each
  eligible structure until no more edges can be oriented. On the other
  hand if the v-structures have been checked previously in either using
  the conservative method or the majority rule, \code{unfVect} contains
  a list of all ambiguous triples in the graph. In this case, only
  structures containing unambiguous triples are oriented until no more
  edges can be oriented. 

  Note that the algorithm works with columns' position of the adjacency
  matrix and not with the names of the variables.

  Note that the order-dependence issues on the 10 orientation rules are
  still present, see Colombo and Maathuis (2013) for more details.
}
\value{
  The outputs are the updated list of sepset and the final adjacency
  matrix M, since some edges may have been deleted. The edge marks in
  the adjacency matrix M are coded in the following way:
  M[i,j]=M[j,i]=0: no edge; M[i,j]=1, M[j,i] != 0: i *-o j; M[i,j]=2,
  M[j,i]!=0: i*->j; M[i,j]=3, M[j,i]!=0: i*-j.   
}
\references{
  D. Colombo and M.H. Maathuis (2013).
  \emph{Order-independent constraint-based causal structure learning},
  (arXiv:1211.3295v2)
  
  D. Colombo, M. H. Maathuis, M. Kalisch, T. S. Richardson
  (2012). Learning high-dimensional directed acyclic graphs with latent
  and selection variables. \emph{Ann. Statist.} \bold{40}, 294-321. 
  
  J. Zhang (2008). On the completeness of orientation rules for causal
  discovery in the presence of latent confounders and selection bias.
  \emph{Artificial Intelligence} \bold{172} 1873-1896.
}
\author{Diego Colombo and Markus Kalisch (\email{kalisch@stat.math.ethz.ch})}

\seealso{\code{\link{rfci}}}

\examples{
##################################################
## Example with hidden variables
## Zhang (2008), Fig. 6, p.1882
##################################################

## draw DAG with latent variables
amat <- t(matrix(c(0,1,0,0,1, 0,0,1,0,0, 0,0,0,1,0, 0,0,0,0,0, 0,0,0,1,0),5,5))
colnames(amat) <- rownames(amat) <- as.character(1:5)
V <- as.character(1:5)
edL <- vector("list",length=5)
names(edL) <- V
edL[[1]] <- list(edges=c(2,4),weights=c(1,1))
edL[[2]] <- list(edges=3,weights=c(1))
edL[[3]] <- list(edges=5,weights=c(1))
edL[[4]] <- list(edges=5,weights=c(1))
g <- new("graphNEL", nodes=V, edgeL=edL,edgemode="directed")
if (require(Rgraphviz)) {
plot(g)
}

## define the latent variable
L <- 1

## compute the true covariance matrix of g
cov.mat <- trueCov(g)
      
## delete rows and columns which belong to L
true.cov <- cov.mat[-L,-L]
      
## transform it in a correlation matrix
true.corr <- cov2cor(true.cov)

## generate 100000 samples of DAG using standard normal error distribution
require(mvtnorm)
n <- 100000
alpha <- 0.01
d.mat <- rmvnorm(n, mean = rep(0,dim(true.corr)[1]), sigma = true.cov) 

## estimate the skeleton of given data
suffStat <- list(C = cor(d.mat), n = n)
indepTest <- gaussCItest
resD <- skeleton(suffStat, indepTest, dim(true.corr)[2], alpha = alpha,
verbose=FALSE)

## estimate all ordered unshielded triple
amat.resD <- as(resD@graph, "matrix")
sepset <- resD@sepset
tmp <- find.unsh.triple(amat.resD,dim(true.corr)[2])
listM <- tmp$unshTripl
vectM <- tmp$unshVect
 
## check and orient v-structures
vstruct.amat <- rfci.vstructures(suffStat, indepTest, dim(true.corr)[2],
alpha = alpha, sepset, amat.resD, listM, vectM, verbose = TRUE)

## estimate the final skeleton and extend it into a PAG using all 10 rules
rules <- rep(TRUE,10)
resP <- udag2apag(vstruct.amat$graph, suffStat, indepTest, alpha = alpha,
vstruct.amat$sepset, rules = rules, verbose = TRUE)
Amat <- resP$graph
colnames(Amat) <- rownames(Amat) <- as.character(2:5)
}
\keyword{multivariate}
\keyword{models}
\keyword{graphs}


