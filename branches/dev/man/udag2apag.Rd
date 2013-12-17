\name{udag2apag}
\alias{udag2apag}
\title{Last step of RFCI algorithm: Transform partially oriented graph into RFCI-PAG}
\description{
  This function performs the last step of the RFCI algorithm: It
  transforms a partially oriented graph in which  
  the v-structures have been oriented into an RFCI Partial Ancestral
  Graph (PAG) (see Colombo et al (2012)).

  While orienting the edges, this function performs some additional
  conditional independence tests in orientation rule 4 to ensure
  correctness of the ancestral relationships.  As a result of these
  additional tests, some additional edges can be deleted.  The result is
  the final adjacency matrix indicating also the edge marks and the
  updated sepsets.
}
\usage{
udag2apag(apag, suffStat, indepTest, alpha, sepset,
          rules = rep(TRUE, 10), unfVect = NULL, verbose = FALSE)
}
\arguments{
  \item{apag}{adjacency matrix of the pcAlgo-object of size p*p where
    the v-structures have already been oriented.}
  \item{suffStat}{sufficient statistics: A \code{\link{list}} containing
    all necessary elements for the conditional independence decisions in the
    function \code{indepTest}.}
  \item{indepTest}{pre-defined function for testing conditional
    independence.  The function is internally called as
    \code{indepTest(x,y,S,suffStat)}, and tests conditional independence
    of \code{x} and \code{y} given \code{S}.  Here, \code{x} and \code{y}
    are variables, and \code{S} is a (possibly empty) set of variables
    (all variables are coded by their column numbers in the adjacency
    matrix).  \code{suffStat} is a list containing all relevant elements
    for the conditional independence decisions.  The return value of 
    \code{indepTest} is the p-value of the test for conditional
    independence.} 
  \item{alpha}{significance level for the individual conditional
    independence tests.}
  \item{sepset}{list of length p; each element of the list
    contains another list of length p.  The element
    \code{sepset[[x]][[y]]} contains the separation set that made the edge
    between \code{x} and \code{y} drop out.  This object is thought to be
    obtained from a \code{pcAlgo}-object.}
  \item{rules}{logical vector of length 10 with \code{TRUE} or \code{FALSE} for 
    each rule, where \code{TRUE} in position i means that rule i (Ri) will be
    applied.  By default, all rules are active.}
  \item{unfVect}{vector containing numbers that encode the ambiguous
    triples (as returned by \code{\link{pc.cons.intern}}).  This is
    needed in the conservative and in the majority rule versions of RFCI.}
  \item{verbose}{logical indicating if detailed output is to be given.}
}
\details{
  The partially oriented graph in which the v-structures have been 
  oriented is transformed into an RFCI-PAG using adapted rules of Zhang
  (2008).  This function is similar to \code{\link{udag2pag}} used to
  orient the skeleton into a PAG in the FCI algorithm.  However, it is
  slightly more complicated because we perform additional conditional
  independence tests when applying rule 4, to ensure correctness of the
  ancestral relationships.  As a result, some additional edges can be
  deleted, see Colombo et al. (2012).  Because of these addiitonal
  tests, we need to give \code{suffStat}, \code{indepTest}, and
  \code{alpha} as inputs. Since edges can be deleted, the input
  adjacency matrix \code{apag} and the input separating sets
  \code{sepset} can change in this algorithm. 

  If \code{unfVect = NULL} (no ambiguous triples), the orientation rules
  are applied to each eligible structure until no more edges can be
  oriented.  On the other hand, hand, if one uses conservative or
  majority rule FCI and ambiguous triples have been found in
  \code{\link{pc.cons.intern}}, \code{unfVect} contains the numbers of
  all ambiguous triples in the graph.  In this case, the orientation
  rules take this information into account.  For example, if a *-> b o-*
  c and <a,b,c> is an unambigous unshielded triple and not a
  v-structure, then we obtain b -* c (otherwise we would create an
  additional v-structure).  On the other hand, if a *-> b o-* c but
  <a,b,c> is an ambiguous unshielded triple, then the circle mark at b
  is not oriented.

  Note that the algorithm works with columns' position of the adjacency
  matrix and not with the names of the variables.
 
  Note that this function does not resolve possible order-dependence in
  the application of the orientation rules, see Colombo and Maathuis (2013).  
}
\value{
  \item{apag}{final adjacency matrix.  The edge marks in
  the adjacency matrix M are coded in the following way:
  M[i,j]=M[j,i]=0: no edge; M[i,j]=1, M[j,i] != 0: i *-o j; M[i,j]=2,
  M[j,i]!=0: i*->j; M[i,j]=3, M[j,i]!=0: i*-j.}

  \item{sepset}{updated list of separating sets} 
}
\references{
  D. Colombo and M.H. Maathuis (2013).
  Order-independent constraint-based causal structure learning.
  arXiv:1211.3295v2.
  
  D. Colombo, M. H. Maathuis, M. Kalisch, T. S. Richardson (2012). 
  Learning high-dimensional directed acyclic graphs with latent
  and selection variables. \emph{Ann. Statist.} \bold{40}, 294--321. 
  
  J. Zhang (2008). On the completeness of orientation rules for causal
  discovery in the presence of latent confounders and selection bias.
  \emph{Artificial Intelligence} \bold{172}, 1873--1896.
}
\author{Diego Colombo and Markus Kalisch (\email{kalisch@stat.math.ethz.ch})}

\seealso{\code{\link{rfci}}, \code{\link{udag2pag}},
  \code{\link{dag2pag}}, \code{\link{udag2pdag}}, 
  \code{\link{udag2pdagSpecial}}, \code{\link{udag2pdagRelaxed}}
}
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
resD <- skeleton(suffStat, indepTest, dim(true.corr)[2],
                 alpha = alpha, verbose=FALSE)

## estimate all ordered unshielded triple
amat.resD <- as(resD@graph, "matrix")
sepset <- resD@sepset
tmp <- find.unsh.triple(amat.resD,dim(true.corr)[2])
listM <- tmp$unshTripl
vectM <- tmp$unshVect
 
## check and orient v-structures
vstrucs <- rfci.vStruc(suffStat, indepTest, p = dim(true.corr)[2],
                       alpha=alpha, sepset=sepset, graph=amat.resD,
                       unshTripl= listM,
                       unshVect = vectM, verbose = TRUE)

## estimate the final skeleton and extend it into a PAG using all 10 rules
rules <- rep(TRUE,10)
resP <- udag2apag(vstrucs$graph, suffStat, indepTest, alpha = alpha,
                  vstrucs$sepset, rules = rules, verbose = TRUE)
Amat <- resP$graph
colnames(Amat) <- rownames(Amat) <- as.character(2:5)
}
\keyword{multivariate}
\keyword{models}
\keyword{graphs}


