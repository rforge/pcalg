\name{udag2pag}
\alias{udag2pag}
\title{Extend the final skeleton into an FCI-PAG}
\description{
This function extends an undirected final skeleton into a Partial
Ancestral Graph (PAG). The final skeleton must have been estimated
either with the function \link{skeleton} (if the search on
Possible-D-SEP has not been performed) or with the
function \link{pdsep}. The result is an adjacency matrix
indicating also the edge marks.
}
\usage{
udag2pag(pag, sepset, rules = rep(TRUE, 10), unfVect = NULL, verbose = FALSE)
}
\arguments{
  \item{pag}{adjacency matrix of the final skeleton of size p*p.}
  \item{sepset}{List of length p; each element of the list
    contains another list of length p. The element
    \code{sepset[[x]][[y]]} contains the separation set that made the edge
    between \code{x} and \code{y} drop out. This object is thought to be
    obtained from a \code{pcAlgo}-object.}
  \item{rules}{array of length 10 containing \code{TRUE} or \code{FALSE}
    for each rule. \code{TRUE} in position i means that rule i (Ri) will
    be used. Per default, all rules are set to TRUE.}
  \item{unfVect}{Vector containing numbers that encode the ambiguous
    triple. This is needed in the conservative and majority rule FCI.}
  \item{verbose}{If \code{TRUE}, detailed output is provided.}
}
\details{
  The skeleton is extended to an FCI-PAG using rules by Zhang (see
  References). Note that the algorithm works with columns' position of the
  adjacency matrix and not with the names of the variables.

  If \code{unfVect = NULL} the v-structures and the ten orientation
  rules are applied to each eligible structure until no more edges can
  be oriented. On the other hand if the v-structures have been checked
  previously either using the conservative method or the majority rule,
  \code{unfVect} contains a list of all ambiguous triples in the
  graph. In this case, only structures containing unambiguous triples
  are oriented until no more edges can be oriented.

  Note that the algorithm works with columns' position of the adjacency
  matrix and not with the names of the variables.

  Note that the order-dependence issues on the 10 orientation rules are
  still present, see Colombo and Maathuis (2013) for more details.
}
\value{
  The output is an adjacency matrix M with edge marks. The edge marks are
  coded in the following way: M[i,j]=M[j,i]=0: no edge; M[i,j]=1, M[j,i]
  != 0: i *-o j; M[i,j]=2, M[j,i]!=0: i*->j; M[i,j]=3, M[j,i]!=0: i*-j.
}
\references{
  D. Colombo and M.H. Maathuis (2013).
  \emph{Order-independent constraint-based causal structure learning},
  (arXiv:1211.3295v2)
  
  J. Zhang (2008). On the completeness of orientation rules for causal
  discovery in the presence of latent confounders and selection bias.
  \emph{Artificial Intelligence} \bold{172} 1873-1896.
}
\author{Diego Colombo and Markus Kalisch (\email{kalisch@stat.math.ethz.ch})}

\seealso{\code{\link{fci}}}

\examples{
##################################################
## Example with hidden variables
## Zhang (2008), Fig. 6, p.1882
##################################################

## draw a DAG with latent variables
## this example is taken from Zhang (2008), Fig. 6, p.1882 (see references)
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

## estimate v-structures conservatively
tmp <- pc.cons.intern(resD, suffStat, indepTest, alpha, verbose = FALSE, version.unf = c(1, 1))
tripleList <- tmp$unfTripl
resD <- tmp$sk

## estimate the final skeleton of given data using Possible-D-Sep
pdsepRes <- pdsep(resD@graph, suffStat, indepTest, dim(true.corr)[2],
resD@sepset, alpha = alpha, m.max = Inf, pMax = resD@pMax, verbose = FALSE)

## extend the skeleton into a PAG using all 10 rules
rules <- rep(TRUE,10)
resP <- udag2pag(pag = pdsepRes$G, pdsepRes$sepset, rules = rules,
verbose = TRUE)
colnames(resP) <- rownames(resP) <- as.character(2:5)
}
\keyword{multivariate}
\keyword{models}
\keyword{graphs}


