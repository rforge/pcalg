\name{dreach}
\alias{dreach}
\title{Compute D-SEP of nodes x and y}
\description{
  This function computes the D-SEP(x,y) set in a maximal ancestral graph
  (MAG). It is intended for internal use only.
}
\usage{
dreach(x, y, amat, verbose = FALSE)
}
\arguments{
\item{x}{Column of node in adjacency matrix, of which D-SEP is
  to be computed}
\item{y}{Column of node in adjacency matrix, of which D-SEP is
  to be computed}
\item{amat}{Adjacency matrix as returned by
  amat[i,j] = 0 iff no edge btw i,j
  amat[i,j] = 2 iff i *-> j
  amat[i,j] = 3 iff i *-- j}
\item{verbose}{Details on output}
}

\value{
  Vector of column positions indicating the nodes in D-SEP(x,y).
}
\references{
  P. Spirtes, C. Glymour and R. Scheines (2000)
  \emph{Causation, Prediction, and Search}, 2nd edition, The MIT Press.
}
\author{Diego Colombo and Markus Kalisch (\email{kalisch@stat.math.ethz.ch})}
\seealso{\code{\link{gen.backdoor}} uses this function.}



