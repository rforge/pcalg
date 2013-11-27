\name{showEdgeList}% FIXME: "too much" in parallel with showAmat.Rd % (notably example)
\alias{showEdgeList}
\title{Show Edge List of pcAlgo object}
\description{
  Show the list of edges (of the graph) of a pcAlgo object; this is
  intended to be an alternative if \pkg{Rgraphviz} does not work.
}
\usage{
showEdgeList(pc.fit, labels = NULL)
}
\arguments{
  \item{pc.fit}{object of class pcAlgo (as returned from
    \code{\link{skeleton}} or \code{\link{pc}}).}
  \item{labels}{optional labels for nodes; by default, the labels from the pcAlgo % 
    object are used.}
}
\value{none; the purpose is in (the side effect of) printing
  the edge list.
}
\author{Markus Kalisch (\email{kalisch@stat.math.ethz.ch})}
\note{
  There is no corresponding function for fciAlgo objects.
}
\seealso{
  \code{\link{iplotPC}} for plotting a pcAlgo object using the package
  \pkg{igraph}; \code{\link{showAmat}} for showing the adjacency matrix of a
  pcAlgo object.
}
\examples{
## Load predefined data
data(gmG)
n <- nrow(gmG$x)
p <- ncol(gmG$x)

## define independence test (partial correlations)
indepTest <- gaussCItest
## define sufficient statistics
suffStat <- list(C = cor(gmG$x), n = n)
## estimate CPDAG
alpha <- 0.01
pc.fit <- pc(suffStat, indepTest, p, alpha, verbose = TRUE)

## Edge list
showEdgeList(pc.fit)

## Adjacency matrix
showAmat(pc.fit)

## Plot using package igraph
if (require(igraph)) {
 ## show estimated CPDAG using igraph
 iplotPC(pc.fit)
}
}
\keyword{graphs}
