\name{iplotPC}
\alias{iplotPC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function for plotting a pcAlgo object using the package igraph}
\description{Function for plotting a pcAlgo object using the package igraph}
\usage{
iplotPC(pc.fit, labels = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pc.fit}{
Object of class pcAlgo (as returned from \code{\link{skeleton}} or \code{\link{pc}})
}
  \item{labels}{
    Optional labels for nodes; by default, the labels from the pcAlgo
    object are used.
}
}
\details{
This function is an alternative to the plot method of pcAlgo objects if
the Rgraphviz package does not work. It extracts the adjacency matrix
and converts it into an igraph object which is then plotted.
}
\value{
  Plot of pcAlgo object.
}
\author{
  Markus Kalisch (\email{kalisch@stat.math.ethz.ch})
}
\note{
Note that this function does not work on fciAlgo objects, because
different edge marks are needed there.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{showEdgeList}} for printing the edge list of a pcAlgo
object; \code{\link{showAmat}} for printing the adjacency matrix of a
pcAlgo object
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
