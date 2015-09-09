\name{displayAmat}
\alias{displayAmat}
\title{Show Adjacency Matrix of pcAlgo or fciAlgo object}
\description{
  Show the adjacency matrix of a \code{"pcAlgo"} object; this is
  intended to be an alternative if the \pkg{Rgraphviz} package does not work.
}
\usage{
displayAmat(obj)
}
\arguments{
  \item{obj}{an \R object of class \code{\linkS4class{pcAlgo}}, as returned
    from \code{\link{skeleton}()} or \code{\link{pc}()} or an object of
    class \code{\linkS4class{fciAlgo}}, as from \code{\link{fci}()}.}
}
\value{
  A \code{\link{list}} with two components:
  \item{amat}{adjacency matrix of type \code{amat.cpdag} if the input
    was an object of class \code{\linkS4class{pcAlgo}}; adjacency matrix
    of type \code{amat.cpdag} if the input was an object of class
    \code{\linkS4class{fciAlgo}}}
  \item{type}{type of adjacency matrix as described in \link{amatType}.}
}
\author{
Markus Kalisch (\email{kalisch@stat.math.ethz.ch})
}
\seealso{
  \code{\link{showEdgeList}} for showing the edge list of a
  \code{\linkS4class{pcAlgo}} object.
  \code{\link{iplotPC}} for plotting a \code{"pcAlgo"} object using
  package \pkg{igraph}.
}
\examples{
##################################################
## pcAlgo object
##################################################
## Load predefined data
data(gmG)
n <- nrow    (gmG8$x)
V <- colnames(gmG8$x)

## define sufficient statistics
suffStat <- list(C = cor(gmG8$x), n = n)
## estimate CPDAG
skel.fit <- skeleton(suffStat, indepTest = gaussCItest,
                     alpha = 0.01, labels = V)
displayAmat(skel.fit)

pc.fit <- pc(suffStat, indepTest = gaussCItest,
             alpha = 0.01, labels = V)
displayAmat(pc.fit)

##################################################
## fciAlgo object
##################################################
set.seed(42)
p <- 7
## generate and draw random DAG :
myDAG <- randomDAG(p, prob = 0.4)

## find PAG using the FCI algorithm
myC <- cov2cor(trueCov(myDAG))
suffStat <- list(C = myC, n = 10^9)
V <- LETTERS[1:p] ## labels of nodes

tmpFCI <- fci(suffStat, indepTest=gaussCItest, labels = V,
              alpha = 0.9999, doPdsep = FALSE)
displayAmat(tmpFCI)
}
\keyword{hplot}