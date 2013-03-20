\name{Score-class}
\docType{class}
\alias{Score-class}
\encoding{UTF-8}
\title{Virtual Class "Score"}
\description{
  This virtual base class represents a score for causal inference; it is used
  in the causal inference functions \code{\link{ges}}, \code{\link{gies}} and
  \code{\link{simy}}. 
}
\section{Fields}{
  The fields of \code{Score} are mainly of interest for users who aim at 
  deriving an own class from this virtual base class, i.e., implementing an own
  score function.
  \describe{
    \item{\code{decomp}:}{Indicates whether the represented score is
      decomposable (cf. details)}
    \item{\code{pp.dat}:}{list representing the preprocessed input data;
      this is typically a statistic which is sufficient for the 
      calculation of the score.}
    \item{\code{.pardag.class}:}{Name of the class of the parametric DAG model
      corresponding to the score.  This must name a class derived from 
      \code{ParDAG}.}
    \item{\code{c.fcn}:}{Only used internally; must remain empty for (user
      specified) classes derived from \code{Score}.} 
  }
}
\section{Methods}{
  Note that since \code{Score} is a virtual class, its methods cannot be called
  directly, but only on derived classes.
  \describe{
    \item{\code{initialize(data = matrix(1, 1, 1), 
      targets = list(integer(0)),
      target.index = rep(as.integer(1), nrow(data)), 
      ...)}}{Constructor; creates a \code{Score} object.  For the meaning of the 
      arguments, see details.}
    \item{\code{local.score(vertex, parents, ...)}}{For decomposable scores,
      this function calculates the local score of a vertex and its parents.
      Must throw an error in derived classes that do not represent a 
      decomposable score.}
    \item{\code{global.score.int(edges, ...)}}{Calculates the global score of
      a DAG, represented as a list of in-edges: for each vertex in the DAG,
      this list contains a vector of parents.}
    \item{\code{global.score(dag, ...)}}{Calculates the global score of a DAG,
      represented as object of a class derived from \code{ParDAG}.}
    \item{\code{local.mle(vertex, parents, ...)}}{Calculates the local MLE of 
      a vertex and its parents.  The result is a vector of parameters whose 
      meaning depends on the model class; it matches the convention used in the
      corresponding causal model (cf. \code{.pardag.class}).}
    \item{\code{global.mle(dag, ...)}}{Calculates the global MLE of a DAG,
      represented by an object of the class specified by \code{.pardag.class}.
      The result is a list of vectors, one per vertex, each in the same format
      as the result vector of \code{local.mle}.}
  }
}
\details{

}
\author{
  Alain Hauser (\email{hauser@stat.math.ethz.ch})
}
\seealso{
  \code{\link{ges}}, \code{\link{gies}}, \code{\link{simy}}, 
  \code{\linkS4class{GaussL0penIntScore}}, 
  \code{\linkS4class{GaussL0penObsScore}}
}
\keyword{classes}
