\name{gmL}
\alias{gmL}
\docType{data}
\title{Latent Variable Graphical Model Data Example}
\description{
This data set contains a matrix containing information on four gaussian
variables and the corresonding DAG model containing four observed and
one latent variable.
}
\usage{data(gmL)}
\format{
  The format is a list of 2 components
  \describe{
    \item{x:}{num [1:10000, 1:4] -0.614 1.304 1.45 -1.121 -0.21 ...}
    \item{g:}{Formal class 'graphNEL' [package "graph"] with 6 slots\cr
      .. ..@ nodes     : chr [1:5] "1" "2" "3" "4" ...\cr
      .. ..@ edgeL     :List of 5\cr
      .......
    }
  }
}
\details{
The data was generated as indicated below. First, a random DAG model was
generated with five nodes; then 10000 samples were drawn from this
model; finally, one of the four variables was declared to be latent and
the corresponding column was deleted from the simulated data set.
}
% \source{
% %%  ~~ reference to a publication or URL from which the data were obtained ~~
% }
% \references{
 %%  ~~ possibly secondary sources and usages ~~
% }
\examples{
data(gmL)
str(gmL, max=3)

## the graph:
gmL$g ; nodes(gmL$g) ; str(edges(gmL$g))
if(require("Rgraphviz"))
  plot(gmL$g, main = "gmL $ g -- latent variable example data")

pairs(gmL $x) # the data
}
\keyword{datasets}
