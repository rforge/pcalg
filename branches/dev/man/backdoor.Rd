\name{backdoor}
\alias{backdoor}
\title{Determine a set that satisfies the Generalized Backdoor Criterion
relative to x and y}
\description{
  This function first checks if the total causal effect of
  one variable (\code{x}) onto another variable (\code{y}) is
  identifiable via the generalized backdoor criterion, and if this is
  the case it explicitly gives a set of variables that satisfies the
  generalized backdoor criterion with respect to \code{x} and \code{y}
  in the given graph.
}
\usage{
backdoor(amat, x, y, type = "pag")
}
\arguments{
  \item{amat}{Adjacency matrix of the given graph specified in \code{type}.}
  \item{x}{Position of variable \code{x} in the adjacency matrix.}
  \item{y}{Position of variable \code{y} in the adjacency matrix.}
  \item{type}{It specifies the type of graph of the adjacency matrix
    \code{amat}. It can be a DAG (type="dag"), a CPDAG (type="cpdag"), a
    MAG (type="mag"), or a PAG (type="pag").}
}
\details{
  This function is a generalization of Pearl's backdoor criterion, see
  Pearl (1993), defined for directed acyclic graphs (DAGs), for single
  interventions and single outcome variable to more general types of
  graphs (CPDAGs, MAGs, and PAGs) that describe Markov equivalence
  classes of DAGs with and without latent variables but without
  selection variables. For more details see Maathuis and Colombo(2013).

  The motivation to find a set W that satisfies the generalized backdoor
  criterion with respect to \code{x} and \code{y}
  in the given graph relies on the result of the generalized backdoor adjustment
  that says: "If a set of variables W satisfies the generalized backdoor
  criterion relative to \code{x} and \code{y} in the given graph, then
  the causal effect of \code{x} on \code{y} is identifiable and is given
  by:" P(Y|do(X = x)) = sum_W P(Y|X,W).P(W). This result allows to write
  post-intervention densities (the one written using Pearl's
  do-calculus) using only observational densities
  estimated from the data.
  
  If the input graph is a DAG (\code{type="dag"}), this function reduces
  to Pearl's backdoor criterion for single interventions and single
  outcome variable, and the parents of \code{x} in the DAG satisfy the
  backdoor criterion unless \code{y} is a parent of \code{x}.

  If the input graph is a CPDAG C (\code{type="cpdag"}), a MAG M
  (\code{type="mag"}), or a PAG P (\code{type="pag"}) (with both M and P
  not allowing selection variables), this function first checks if the
  total causal effect of \code{x} on \code{y} is identifiable via the
  generalized backdoor criterion (see Maathuis and Colombo, 2013). If
  the effect is not identifiable in this way, the output is
  NA. Otherwise, an explicit set W that satisfies the generalized
  backdoor criterion with respect to \code{x} and \code{y} in the given
  graph is found. 

  At this moment this function is not able to work with an RFCI-PAG. 
  
  It is important to note that there can be pair of nodes \code{x} and
  \code{y} for which there is no set W that satisfies the generalized
  backdoor criterion, but the total causal effect might be identifiable
  via some other technique.
}
\value{
  Either NA if the total causal effect is not identifiable via the
  generalized backdoor criterion, or a set if the effect is identifiable
  via the generalized backdoor criterion. Note that if the set W is
  equal to the empty set, the output is NULL.
}
\references{
  M.H. Maathuis and D. Colombo (2013). A generalized backdoor
  criterion. arXiv preprint arXiv:1307.5636.

  J. Pearl (1993). Comment: Graphical models, causality and intervention.
  \emph{Statistical Science} \bold{8}, 266--269.
}
\author{Diego Colombo and Markus Kalisch (\email{kalisch@stat.math.ethz.ch})}
\seealso{\code{\link{pc}} for estimating a CPDAG, \code{\link{dag2pag}}
  and \code{\link{fci}} for estimating a PAG, and \code{\link{pag2mag}}
  for estimating a MAG.
}
\examples{
#####################################################################
##DAG
#####################################################################
## Simulate the true DAG
set.seed(123)
p <- 7
myDAG <- randomDAG(p, prob = 0.2) ## true DAG

## Extract the adjacency matrix of the true DAG
true.amat <- as(myDAG, "matrix")
true.amat[which(true.amat!=0)] <- 1

## Compute the effect using the generalized backdoor criterion
%% library(igraph)
backdoor(true.amat, 5, 7, type="dag")

#####################################################################
##CPDAG
#####################################################################
##################################################
## Example not identifiable
## Maathuis and Colombo (2013), Fig. 3, p.14
##################################################

## create the graph
p <- 5
amat <- t(matrix(c(0,0,1,1,1, 0,0,1,1,1, 0,0,0,1,0, 0,0,0,0,1, 0,0,0,0,0),5,5))
colnames(amat) <- rownames(amat) <- as.character(1:5)
V <- as.character(1:5)
edL <- vector("list",length=5)
names(edL) <- V
edL[[1]] <- list(edges=c(3,4,5),weights=c(1,1,1))
edL[[2]] <- list(edges=c(3,4,5),weights=c(1,1,1))
edL[[3]] <- list(edges=4,weights=c(1))
edL[[4]] <- list(edges=5,weights=c(1))
g <- new("graphNEL", nodes=V, edgeL=edL, edgemode="directed")

## estimate the true CPDAG
myCPDAG <- dag2cpdag(g)
## Extract the adjacency matrix of the true CPDAG
true.amat <- as(myCPDAG, "matrix")
true.amat[which(true.amat!=0)] <- 1

## The effect is not identifiable, in fact:
%% library(igraph)
backdoor(true.amat, 3, 5, type="cpdag")


##################################################
## Example identifiable
## Maathuis and Colombo (2013), Fig. 4, p.15
##################################################

## create the graph
p <- 6
amat <- rbind(c(0,0,1,1,0,1), c(0,0,1,1,0,1), c(0,0,0,0,1,0),
              c(0,0,0,0,1,1), c(0,0,0,0,0,0), c(0,0,0,0,0,0))
colnames(amat) <- rownames(amat) <- as.character(1:6)
V <- as.character(1:6)
edL <- vector("list",length=6)
names(edL) <- V
edL[[1]] <- list(edges=c(3,4,6),weights=c(1,1,1))
edL[[2]] <- list(edges=c(3,4,6),weights=c(1,1,1))
edL[[3]] <- list(edges=5,weights=c(1))
edL[[4]] <- list(edges=c(5,6),weights=c(1,1))
g <- new("graphNEL", nodes=V, edgeL=edL, edgemode="directed")

## estimate the true CPDAG
myCPDAG <- dag2cpdag(g)
## Extract the adjacency matrix of the true CPDAG
true.amat <- as(myCPDAG, "matrix")
true.amat[which(true.amat!=0)] <- 1

## The effect is identifiable and
%% library(igraph)
backdoor(true.amat, 6, 3, type="cpdag")

##################################################################
##PAG
##################################################################
##################################################
## Example identifiable
## Maathuis and Colombo (2013), Fig. 7, p.17
##################################################

## create the graph
p <- 7
amat <- t(matrix(c(0,0,1,1,0,0,0, 0,0,1,1,0,0,0, 0,0,0,1,0,1,0, 
                   0,0,0,0,0,0,1, 0,0,0,0,0,1,1, 0,0,0,0,0,0,0,
                   0,0,0,0,0,0,0),  7, 7))
colnames(amat) <- rownames(amat) <- as.character(1:7)
V <- as.character(1:7)
edL <- vector("list",length=7)
names(edL) <- V
edL[[1]] <- list(edges=c(3,4),weights=c(1,1))
edL[[2]] <- list(edges=c(3,4),weights=c(1,1))
edL[[3]] <- list(edges=c(4,6),weights=c(1,1))
edL[[4]] <- list(edges=7,weights=c(1))
edL[[5]] <- list(edges=c(6,7),weights=c(1,1))
g <- new("graphNEL", nodes=V, edgeL=edL, edgemode="directed")
L <- 5

## compute the true covariance matrix of g
cov.mat <- trueCov(g)

## transform covariance matrix into a correlation matrix
true.corr <- cov2cor(cov.mat)
suffStat <- list(C=true.corr, n=10^9)
indepTest <- gaussCItest
%% library(RBGL)

## estimate the true PAG
true.pag <- dag2pag(suffStat, indepTest, g, L, alpha = 0.9999)

## The effect is identifiable and
%% library(igraph)
backdoor(true.pag@amat, 4, 6, type="pag")
}
\keyword{multivariate}
\keyword{models}
\keyword{graphs}
