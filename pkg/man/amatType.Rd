\name{amatType}
\alias{amatType}
\alias{amat.cpdag}
\alias{amat.pag}
\title{Types of Adjacency Matrices}
\description{
  Two types of adjacency matrices are used in package \pkg{pcalg}: Type
  \code{amat.cpdag} for DAGs and CPDAGs and type \code{amat.pag} for
  MAGs and PAGs. The required type of adjacency matrix is documented
  in the help files of the respective functions or classes. If in some
  functions more detailed information on the graph type is needed
  (i.e. DAG or CPDAG; MAG or PAG) this information will be passed in a
  separate argument (see e.g. \code{\link{gac}} and the examples below).
}
\details{
  Adjacency matrices are square integer matrices with zeros on the
  diagonal. They can have row- and columnnames; however, most functions
  will work on the (integer) column positions in the adjacency matrix.
  \cr
  Coding for type \code{amat.cpdag}:
  Coding 0: No edge or tail
  Coding 1: Arrowhead
  Note that the edgemark-code refers to the \emph{row}index (as opposed adjacency
  matrices of type mag or pag).  
  E.g.:
  \code{amat[a,b] = 0} and \code{amat[b,a] = 1} implies a --> b.
  \code{amat[a,b] = 1} and \code{amat[b,a] = 0} implies a <-- b.
  \code{amat[a,b] = 0} and \code{amat[b,a] = 0} implies a     b.
  \code{amat[a,b] = 1} and \code{amat[b,a] = 1} implies a --- b.
  
  \cr
  Coding for type \code{amat.pag}:
  Coding 0: No edge
  Coding 1: Circle
  Coding 2: Arrowhead
  Coding 3: Tail
  Note that the edgemark-code refers to the
  \emph{colum}index (as opposed adjacency matrices of type dag or cpdag)
  E.g.:
  \code{amat[a,b] = 2} and \code{amat[b,a] = 3} implies a --> b.
  \code{amat[a,b] = 3} and \code{amat[b,a] = 2} implies a <-- b.
  \code{amat[a,b] = 2} and \code{amat[b,a] = 2} implies a <-> b.
  \code{amat[a,b] = 1} and \code{amat[b,a] = 3} implies a --o b.
  \code{amat[a,b] = 0} and \code{amat[b,a] = 0} implies a     b.
}

\seealso{E.g. \code{\link{gac}} for a function which takes an
  adjacency matrix as input; \code{\link{displayAmat}} for a function which
  returns an 
  adjacency matrix as output; \code{\linkS4class{fciAlgo}} for a class
  which has an  adjacency matrix in one slot.}
\examples{
##################################################
## Function gac() takes an adjecency matrix of
## any kind as input. In addition to that, the
## precise type of graph (DAG/CPDAG/MAG/PAG) needs
## to be passed as a different argument
##################################################
## Adjacency matrix of type \code{amat.cpdag}
m1 <- matrix(c(0,1,0,1,0,0, 0,0,1,0,1,0, 0,0,0,0,0,1,
                  0,0,0,0,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0),6,6)
## more detailed information on the graph type needed by gac()
type <- "dag"
gac(m1, x=1,y=3, z=NULL, type = type)

## Adjacency matrix of type \code{amat.cpdag}
m2 <- matrix(c(0,1,1,0,0,0, 1,0,1,1,1,0, 0,0,0,0,0,1,
                  0,1,1,0,1,1, 0,1,0,1,0,1, 0,0,0,0,0,0), 6,6)
## more detailed information on the graph type needed by gac()
type <- "cpdag"
gac(m2, x=3, y=6, z=c(2,4), type = type)

## Adjacency matrix of type \code{amat.pag}
m3 <- matrix(c(0,2,0,0, 3,0,3,3, 0,2,0,3, 0,2,2,0), 4,4)
## more detailed information on the graph type needed by gac()
type <- "mag"
gac(m3, x=2, y=4, z=NULL,   type = type)
## more detailed information on the graph type needed by gac()
type <- "pag"
gac(m3, x=2, y=4, z=NULL,   type = type)

##################################################
## Function displayAmat() returns an adjacency matrix
## and the type of the adjacency matrix.
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
resSkel <- displayAmat(skel.fit)

##################################################
## Function fci() returns an adjacency matrix
## of type amat.pag as one slot.
##################################################
set.seed(42)
p <- 7
## generate and draw random DAG :
myDAG <- randomDAG(p, prob = 0.4)

## find skeleton and PAG using the FCI algorithm
suffStat <- list(C = cov2cor(trueCov(myDAG)), n = 10^9)
res <- fci(suffStat, indepTest=gaussCItest,
           alpha = 0.9999, p=p, doPdsep = FALSE)
str(res)
}
\keyword{ graphs }
