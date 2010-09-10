\name{gmL}
\alias{gmL}
\docType{data}
\title{Latent Variable Graphical Model Data Example}
\description{
  ......
}
\usage{data(gmL)}
\format{
  The format is:
%% become shorter, use \item -- FIXME
List of 2
 $ x: num [1:10000, 1:4] -0.614 1.304 1.45 -1.121 -0.21 ...

 $ g:Formal class 'graphNEL' [package "graph"] with 6 slots
  .. ..@ nodes     : chr [1:5] "1" "2" "3" "4" ...
  .. ..@ edgeL     :List of 5
  .. ..@ edgeData  :Formal class 'attrData' [package "graph"] with 2 slots
  .. ..@ nodeData  :Formal class 'attrData' [package "graph"] with 2 slots
  .. ..@ renderInfo:Formal class 'renderInfo' [package "graph"] with 4 slots
  .. ..@ graphData :List of 1
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
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
