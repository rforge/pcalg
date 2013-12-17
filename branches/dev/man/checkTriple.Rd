\name{checkTriple}
\alias{checkTriple}
\title{Check Consistency of Conditional Independence for a Triple of Nodes}
\description{
  For each subset of \code{nbrsA} and \code{nbrsC} where \code{a} and
  \code{c} are conditionally independent, it is checked if \code{b} is in the
  conditioning set.  
}
\usage{
checkTriple(a, b, c, nbrsA, nbrsC,
            sepsetA, sepsetC,
            suffStat, indepTest, alpha, version.unf = c(NA, NA),
            maj.rule = FALSE, verbose = FALSE)
}
\arguments{
  \item{a, b, c}{integer positions in adjacency matrix for nodes
    \eqn{a}, \eqn{b}, and \eqn{c}, respectively.}
  \item{nbrsA, nbrsC}{neighbors of \eqn{a} and \eqn{c}, respectively.}
  \item{sepsetA}{vector containing \eqn{Sepset(a,c)}.}
  \item{sepsetC}{vector containing \eqn{Sepset(c,a)}.}
  \item{suffStat}{a \code{\link{list}} of sufficient statistics for
    independent tests; see, e.g., \code{\link{pc}}.}
  \item{indepTest}{a \code{\link{function}} for the independence test,
    see, e.g., \code{\link{pc}}.}
  \item{alpha}{significance level of test.}
  \item{version.unf}{
    Vector of length two:
    \code{version.unf[1]}: 1 - check for all separating subsets of
    \code{nbrsA} and \code{nbrsC} if \code{b} is in that set, 2 - it
    also checks if there at all exists 
    any sepset which is a subset of the neighbours (there might be none,
    although \code{b} is in the sepset, which indicates an ambiguous situation);
    version.unf[2]: 1 - do not consider the initial sepsets
    \code{sepsetA} and \code{sepsetC} (same as Tetrad), 2 - consider if
    \code{b} is in \code{sepsetA} or \code{sepsetC}.}
  \item{maj.rule}{logical, indicating that the following majority rule
    is applied: if \code{b} is in less than 
    50\% of the checked sepsets, we say that \code{b} is in \bold{no} sepset.  If
    \code{b} is in more than 50\% of the checked sepsets, we say that
    \code{b} is in \bold{all} sepsets.  If \code{b} is in exactly 50\% of the
    checked sepsets, the triple is considered \sQuote{ambiguous}.}
  \item{verbose}{logical asking for detailed output of intermediate steps.}
}
\details{
This function is used in the conservative versions of structure learning
  algorithms. 
}
\value{
  \item{decision}{Decision on possibly ambiguous triple with codes
    \describe{
      \item{1}{\code{b} is in NO sepset (make v-structure);}
      \item{2}{\code{b} is in ALL sepsets (make no v-structure);}
      \item{3}{\code{b} is in SOME but not all sepsets (ambiguous triple)}
  }}
  \item{vers}{version (1 or 2) of the ambiguous triple
    (1=normal ambiguous triple that is \code{b} is in some sepsets;
    2=triple coming from \code{version.unf[1]==2}, that is, \code{a} and
    \code{c} are indep given the initial sepset but there doesn't exist a
    subset of the neighbours that d-separates them.)}
  \item{sepsetA}{Updated version of \code{sepsetA}}
  \item{sepsetC}{Updated version of \code{sepsetC}}
}
\references{
  D. Colombo and M.H. Maathuis (2013).
  Order-independent constraint-based causal structure learning,
  (arXiv:1211.3295v2).
}
\author{
  Markus Kalisch (\email{kalisch@stat.math.ethz.ch}) and Diego Colombo.
}

\examples{
##################################################
## Using Gaussian Data
##################################################
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
if (require(Rgraphviz)) {
  ## show estimated CPDAG
  par(mfrow=c(1,2))
  plot(pc.fit, main = "Estimated CPDAG")
  plot(gmG$g, main = "True DAG")
}

a <- 6
b <- 1
c <- 8
nbrsA <- c(1,5,7)
nbrsC <- c(1,5)
sepsetA <- pc.fit@sepset[[a]][[c]]
sepsetC <- pc.fit@sepset[[c]][[a]]

checkTriple(a, b, c, nbrsA, nbrsC, sepsetA, sepsetC, suffStat,
indepTest, alpha, verbose = TRUE, version.unf = c(2,2), maj.rule = FALSE)

}
\keyword{manip}
