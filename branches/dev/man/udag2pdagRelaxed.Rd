\name{udag2pdagRelaxed}
\alias{udag2pdagRelaxed}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Extend a pcAlgo-object containing a skeleton to a PDAG}
\description{
This function extends a pcAlgo-object containing a skeleton and
corresponding conditional independence information to its corresponding
CPDAG. The result is a pcAlgo-object as well. In the oracle version, the
algorithm returns the true CPDAG. There is no check whether the result
is a valid CPDAG, and in the sample version it may be invalid.
}
\usage{
udag2pdagRelaxed(gInput, verbose, unfVect=NULL, solve.confl=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gInput}{pcAlgo-object containing skeleton and
    cond. ind. information}
  \item{verbose}{0: No output; 1: Details}
  \item{unfVect}{Vector containing numbers that encode the ambiguous
    triples (as returned by \code{\link{pc.cons.intern}}). This is
    needed in the conservative and majority rule PC.} 
  \item{solve.confl}{If \code{TRUE}, the orientation of the v-structures and the
    orientation rules work with lists for candidate sets and allow
    bi-directed edges to resolve conflicting edge orientations. Note that
    therefore the resulting object might not be a PDAG because
    bi-directed edges might be present.}
}
\details{
  The skeleton is extended to a PDAG using rules by Meek (see
  References).

  If \code{unfVect = NULL} the three orientation rules are applied to each
  eligible structure until no more edges can be oriented. On the other
  hand if the v-structures have been checked previously in
  \code{\link{pc.cons.intern}} either using the conservative method or
  the majority rule, \code{unfVect} contains a list of all ambiguous
  triples in the graph. In this case, only structures containing
  unambiguous triples are oriented using the three orientation rules
  until no more edges can be oriented. 

  If \code{solve.confl = TRUE} both the v-structures and the orientation
  rules work with lists for the candidate edges and allow bi-directed
  edges if a conflicting orientation is found on the edge. We denote
  bi-directed edges, for example between two variables i and j, in the
  adjacency matrix M of the graph as M[i,j]=2 and M[j,i]=2. In this case
  the output cannot be interpreted causally, because bi-directed edges
  have no causal meaning in a CPDAG, see Colombo and Maathuis (2013).
  
  There is no test whether the result is a valid PDAG. 
}
\value{
  \item{pcObj}{Oriented pc-Object}
}
\references{
  C. Meek (1995). Causal inference and causal explanation with
  background knowledge. In \emph{Proceedings of the Eleventh Conference
    on Uncertainty in Artificial Intelligence (UAI-95)},
  pp. 403–411. Morgan Kaufmann Publishers, Inc.

  P. Spirtes, C. Glymour and R. Scheines (2000)
  \emph{Causation, Prediction, and Search}, 2nd edition, The MIT Press.

  D. Colombo and M.H. Maathuis (2013).
  \emph{Order-independent constraint-based causal structure learning},
  (arXiv:1211.3295v1).
}
\author{Markus Kalisch (\email{kalisch@stat.math.ethz.ch}) and Diego Colombo}
\seealso{\code{\link{pdag2dag}}, \code{\link{dag2cpdag}},
  \code{\link{udag2pdag}}, \code{\link{udag2pdagRelaxed}}, \code{\link{udag2pdagSpecial}}}
\examples{
## simulate data
set.seed(123)
p <- 10
myDAG <- randomDAG(p, prob = 0.2)
trueCPDAG <- dag2cpdag(myDAG)
n <- 1000
library(RBGL)
d.mat <- rmvDAG(n, myDAG, errDist = "normal")

## define independence test (partial correlations)
indepTest <- gaussCItest 
## define sufficient statistics
suffStat <- list(C = cor(d.mat), n = n)
## estimate skeleton
resU <- skeleton(suffStat, indepTest, p, 0.05)

## orient edges using three different methods
resD1 <- udag2pdagRelaxed(resU,verbose=0)
resD2 <- udag2pdagSpecial(resU,verbose=0,n.max=100)
resD3 <- udag2pdag(resU,verbose=0)

}
\keyword{multivariate}
\keyword{models}
\keyword{graphs}


