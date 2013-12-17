\name{find.unsh.triple}
\alias{find.unsh.triple}
\title{
Find all unshielded triples in an undirected graph.
}
\description{
Find the ordered (<x,y,z> with x<z) list of all the triples in the graph.
}
\usage{
find.unsh.triple(g, p)
}
\arguments{
  \item{g}{Adjacency matrix of skeleton
}
  \item{p}{Number of nodes of the graph
}
}
\details{
A triple of nodes \code{x}, \code{y} and \code{z} is unshielded, if: (i)
\code{x} and \code{y} are connected; (ii) \code{y} and \code{z} are
connected; (iii)  \code{x} and \code{z} are NOT connected.
}
\value{
  \item{unshTripl}{Matrix with 3 rows containing in each column
                   an unshielded triple}
  \item{unshVect}{Vector containing the unique number  for
    each column in unshTripl (for internal use only)}
}
\author{
 Diego Colombo and Markus Kalisch (\email{kalisch@stat.math.ethz.ch})
}

\examples{
data(gmG)
if (require(Rgraphviz)) {
  ## show graph
  plot(gmG$g, main = "True DAG")
}

## prepare skeleton use in example
g <- wgtMatrix(gmG$g) ## compute weight matrix
g[which(g != 0)] <- 1 ## compute adjacency matrix
skel <- g + t(g) ## adjacency matrix of skeleton
p <- ncol(skel) ## number of nodes in graph

## estimate unshielded triples
find.unsh.triple(skel, p)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}
