\name{gmL}
\alias{gmL}
\docType{data}
\title{Latent Variable Graphical Model Data Example}
\description{
  ......
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
%%  ~~ If necessary, more details than the __description__ above ~~
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
