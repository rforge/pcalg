\name{pag2mag}
\alias{pag2mag}
\encoding{UTF8}
\title{Transform a PAG into a valid MAG in the Markov equivalence class
  represented by the PAG} 
\description{
  Transform a Partial Ancestral Graph (PAG) into a valid Maximal
  Ancestral Graph (MAG) that belongs to the equivalence class
  represented by the given PAG.
}
\usage{
pag2mag(amat.pag, x)
}
\arguments{
  \item{amat.pag}{Adjacency matrix corresponding to the PAG.}
  \item{x}{Node in the PAG for which no additional edges are oriented
    into it.}
}
\details{
  This function converts a PAG (adjacency matrix) to a valid MAG
  (adjacency matrix) that belongs to the equivalence class
  represented by the given PAG. Note that we assume that there are no
  selection variables, meaning that the edges in the PAG can be of the
  following types: ->, <->, o->, and o-o. In a first step, it uses the
  Arrowhead Augmentation of Zhang (2006), i.e. any o-> edge is oriented
  into ->. Afterwards, it orients any chordal component into a valid DAG
  without orienting any additional edges into x. 
}
\value{
  The output is an adjacency matrix M of a valid MAG with edge
  marks. The edge marks are coded in the following way: M[i,j]=M[j,i]=0:
  no edge; M[i,j]=2, M[j,i]=3: i -> j; M[i,j]=2, M[j,i]=2: i <-> j.
}
\references{
  Zhang, J. (2006). Causal Inference and Reasoning in Causally
  Insufficient Systems. Ph. D. thesis, Carnegie Mellon University.
}
\author{Diego Colombo and Markus Kalisch (\email{kalisch@stat.math.ethz.ch}).}
\seealso{\code{\link{fci}}, \code{\link{dag2pag}}}
\examples{
## create the graph
set.seed(78)
p <- 12
g <- randomDAG(p, prob = 0.4)

## define nodes 2 and 6 to be latent variables
L <- c(2,6)

## compute the true covariance matrix of g
cov.mat <- trueCov(g)
      
## transform covariance matrix into a correlation matrix
true.corr <- cov2cor(cov.mat)

## find PAG
## as dependence "oracle", we use the true correlation matrix in the
## function gaussCItest with a large "virtual sample size" and a large
## alpha
suffStat <- list(C=true.corr, n=10^9)
indepTest <- gaussCItest
library(RBGL)
true.pag <- dag2pag(suffStat, indepTest, g, L, alpha = 0.9999)

## find a valid MAG such that no additional edges are directed into
library(igraph)
valid.mag <- pag2mag(true.pag@amat, 4)
}
\keyword{multivariate}
\keyword{models}
\keyword{graphs}

