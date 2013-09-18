\name{dag2cpdag}
\alias{dag2cpdag}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Convert a DAG to a CPDAG}
\description{
Convert a DAG to a Completed Partially Directed Acyclic Graph (CPDAG).
}
\usage{
dag2cpdag(g)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{g}{DAG (graph object)}
}
\details{
  This function converts a DAG (graph object) into its corresponding
  (unique) CPDAG (graph object) as follows. Because every DAG in the
  Markov equivalence class described by a CPDAG shares the same skeleton
  and the same v-structures, this function takes the skeleton and the
  v-structures of the given DAG \code{g}. Afterwards it simply uses the
  3 orientation rules of the PC algorithm (see references) to orient as
  many of the remaining undirected edges as possible. 
}
\value{
  A graph object containing the CPDAG.
}
\references{
  C. Meek (1995). Causal inference and causal explanation with
  background knowledge. In \emph{Proceedings of the Eleventh Conference on
  Uncertainty in Artificial Intelligence (UAI-95)}, pp. 403–411. Morgan
  Kaufmann Publishers, Inc.

  P. Spirtes, C. Glymour and R. Scheines (2000)
  \emph{Causation, Prediction, and Search}, 2nd edition, The MIT Press.
}
\author{Markus Kalisch (\email{kalisch@stat.math.ethz.ch}) and Diego Colombo}
\seealso{\code{\link{udag2pdag}}, \code{\link{pdag2dag}}}
\examples{
p <- 10 # number of random variables
s <- 0.4 # sparseness of the graph

## generate random data
set.seed(42)
g <- randomDAG(p,s) # generate a random DAG

res <- dag2cpdag(g)
}
\keyword{multivariate}
\keyword{models}
\keyword{graphs}

