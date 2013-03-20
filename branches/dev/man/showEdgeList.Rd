\name{showEdgeList}
\alias{showEdgeList}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Show edge list of pcAlgo object
}
\description{
  Show edge list of pcAlgo object; this is intended to be an alternative
  if Rgraphviz does not work.
}
\usage{
showEdgeList(pc.fit, labels = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pc.fit}{
Object of class pcAlgo (as returned from \code{\link{skeleton}} or \code{\link{p
c}})
}
  \item{labels}{
    Optional labels for nodes; by default, the labels from the pcAlgo
    object are used.
}
}
\value{
Print out of edge list.
}
\author{
Markus Kalisch (\email{kalisch@stat.math.ethz.ch})
}
\note{
There is no correspondin function for fciAlgo objects.
}
\seealso{
\code{\link{iplotPC}} for plotting a pcAlgo object using the package
igraph; \code{\link{showAmat}} for showing the adjacency matrix of a
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{graphs}
