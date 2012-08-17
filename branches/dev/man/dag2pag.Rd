\name{dag2pag}
\alias{dag2pag}

\title{Convert a DAG with latent variables to a PAG }
\description{
Convert a DAG with latent variables to a Partially Ancestral Graph (PAG).
}
\usage{
dag2pag(suffStat, indepTest, graph, ancList, L, alpha, rules = rep(TRUE,10), verbose = FALSE)
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
  \item{g}{DAG with \code{p} nodes (graph object).}
  \item{ancList}{List of length \code{p} containing the ancestors of
  each node in the graph \code{g}. This must be precomputed using the
  function ancTS(g).}
  \item{L}{Array containing the labels of the nodes in the graph
  \code{g} corresponding to the latent variables.}
  \item{alpha}{Significance level for the individual conditional
    independence tests.}
  \item{rules}{Logical vector of length 10 indicating which rules
    should be used when directing edges. The order of the rules is taken
    from Zhang (2009).}
  \item{verbose}{If \code{TRUE}, detailed output is provided.}
}
\details{
  This function converts a DAG (graph object) with latent variables to
  its corresponding (unique) PAG (fciAlgo object), using the ancestors
  information entailed in the true DAG. The output of this function is
  exactly the same as the one using fci(suffStat, indepTest, p, alpha,
  rules = rep(TRUE, 10), doPdsep = TRUE), using the true correlation
  matrix in the function \code{gaussCItest} with a large "virtual sample
  size" and a large alpha, but much more faster.
}
\value{An object of \code{\link{class}} \code{fciAlgo} (see
 \code{\linkS4class{fciAlgo}}) containing the estimated graph
 (in the form of an adjacency matrix with various possible edge marks),
 the conditioning sets that lead to edge removals (sepset) and several other
 parameters.
}
\references{
  Richardson, T. and Spirtes, P. (2002). Ancestral graph Markov
  models. \emph{Ann. Statist.} \bold{30}, 962–1030, Theorem 4.2., page 983.
}
\author{Diego Colombo and Markus Kalisch (\email{kalisch@stat.math.ethz.ch}).}
\seealso{\code{\link{fci}}, \code{\link{pc}}}
\examples{
## create the graph
set.seed(78)
p <- 8
g <- randomDAG(p, prob = 0.4)

## compute the ancestors
library(RBGL)
anc.list <- ancTS(g)

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
true.pag <- dag2pag(suffStat, indepTest, g, anc.list, L, alpha=0.9999)

##find PAG using fci-function
## compute the true covariance matrix of g1
cov.mat1 <- trueCov(g)

## delete rows and columns belonging to latent variable L
true.cov1 <- cov.mat1[-L,-L]

## transform covariance matrix into a correlation matrix
true.corr1 <- cov2cor(true.cov1)

## find PAG with FCI algorithm
## as dependence "oracle", we use the true correlation matrix in the
## function gaussCItest with a large "virtual sample size" and a large
## alpha
suffStat1 <- list(C = true.corr1, n = 10^9)
indepTest1 <- gaussCItest
true.pag1 <- fci(suffStat1, indepTest1, dim(true.corr1)[1], alpha = 0.9999)

##check if the outputs are equal
any(true.pag@amat != true.pag1@amat)
}
\keyword{multivariate}
\keyword{models}
\keyword{graphs}

