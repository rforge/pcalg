\name{showAmat}
\alias{showAmat}
\title{Show Adjacency Matrix of pcAlgo object}
\description{
  Show the  adjacency matrix of a \code{"pcAlgo"} object; this is
  intended to be an alternative if the \pkg{Rgraphviz} package does not work.
}
\usage{
showAmat(pc.fit)
}
\arguments{
  \item{pc.fit}{
    Object of class \code{"pcAlgo"} (as returned from
    \code{\link{skeleton}} or \code{\link{pc}}).}
}
\value{
  The adjacency matrix.
}
\author{
Markus Kalisch (\email{kalisch@stat.math.ethz.ch})
}
\note{
  For \code{"fciAlgo"} objects, the show method produces a similar result.
}

\seealso{
  \code{\link{iplotPC}} for plotting a \code{"pcAlgo"} object using the package
  \pkg{igraph}; \code{\link{showEdgeList}} for showing the edge list of a
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

## Plot: show estimated CPDAG using igraph
iplotPC(pc.fit)
}
\keyword{hplot}
