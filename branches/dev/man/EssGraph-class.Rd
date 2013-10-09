\name{EssGraph-class}
\docType{class}
\alias{EssGraph-class}
\alias{plot,EssGraph,ANY-method}

\title{Class \code{"EssGraph"}}
\description{
  This class represents an (observentional or interventional) essential graph.
}
\section{Extends}{

All reference classes extend and inherit methods from \code{"\linkS4class{envRefClass}"}.

}
\section{Methods}{
  \describe{
    \item{plot}{\code{signature(x = "EssGraph", y = "ANY")}: ... }
	 }
}
\author{
  Alain Hauser (\email{alain.hauser@biology.unibe.ch})
}
\seealso{
  \code{\linkS4class{ParDAG}}
}
\examples{
showClass("EssGraph")
}
\keyword{classes}
\section{Fields}{
  \describe{
    \item{\code{.nodes}:}{Object of class \code{vector} ~~ }
    \item{\code{.in.edges}:}{Object of class \code{list} ~~ }
    \item{\code{targets}:}{Object of class \code{list} ~~ }
    \item{\code{score}:}{Object of class \code{Score} ~~ }
  }
}
\section{Class-Based Methods}{
  \describe{
    \item{\code{greedy.step(direction)}:}{ ~~ }
    \item{\code{repr()}:}{ ~~ }
    \item{\code{causal.inf.options(caching, turning, maxdegree, maxsteps, childrenonly, DEBUG.LEVEL)}:}{ ~~ }
    \item{\code{mle.fit(dag)}:}{ ~~ }
    \item{\code{opt.target(max.size)}:}{ ~~ }
    \item{\code{caus.inf(algorithm, ...)}:}{ ~~ }
    \item{\code{node.count()}:}{ ~~ }
    \item{\code{initialize(nodes, in.edges, ...)}:}{ ~~ }
    \item{\code{edge.count()}:}{ ~~ }
    \item{\code{greedy.search(direction)}:}{ ~~ }
    \item{\code{silander(...)}:}{ ~~ }
    \item{\code{gies(...)}:}{ ~~ }
    \item{\code{gds(...)}:}{ ~~ }
  }
}
