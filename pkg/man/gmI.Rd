\name{gmI}
\alias{gmI}
\docType{data}
\title{Graphical Model IDA Data Example}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(gmI)}
\format{
  The format is:
List of 2
 $ x: num [1:10000, 1:7] -1.727 1.69 0.504 2.528 0.549 ...
 $ g:Formal class 'graphNEL' [package "graph"] with 6 slots
  .. ..@ nodes     : chr [1:7] "1" "2" "3" "4" ...
  .. ..@ edgeL     :List of 7
  .. .. ..$ 1:List of 1
  .. .. .. ..$ edges: int 6
  .. .. ..$ 2:List of 1
  .. .. .. ..$ edges: int [1:2] 7 3
  .. .. ..$ 3:List of 1
  .. .. .. ..$ edges: int 5
  .. .. ..$ 4:List of 1
  .. .. .. ..$ edges: int(0)
  .. .. ..$ 5:List of 1
  .. .. .. ..$ edges: int 7
  .. .. ..$ 6:List of 1
  .. .. .. ..$ edges: num 7
  .. .. ..$ 7:List of 1
  .. .. .. ..$ edges: int(0)
  .. ..@ edgeData  :Formal class 'attrData' [package "graph"] with 2 slots
  .. .. .. ..@ data    :List of 6
  .. .. .. .. ..$ 1|6:List of 1
  .. .. .. .. .. ..$ weight: num 0.468
  .. .. .. .. ..$ 2|7:List of 1
  .. .. .. .. .. ..$ weight: num 0.575
  .. .. .. .. ..$ 2|3:List of 1
  .. .. .. .. .. ..$ weight: num 0.903
  .. .. .. .. ..$ 3|5:List of 1
  .. .. .. .. .. ..$ weight: num 0.961
  .. .. .. .. ..$ 5|7:List of 1
  .. .. .. .. .. ..$ weight: num 0.193
  .. .. .. .. ..$ 6|7:List of 1
  .. .. .. .. .. ..$ weight: num 0.321
  .. .. .. ..@ defaults:List of 1
  .. .. .. .. ..$ weight: num 1
  .. ..@ nodeData  :Formal class 'attrData' [package "graph"] with 2 slots
  .. .. .. ..@ data    : list()
  .. .. .. ..@ defaults: list()
  .. ..@ renderInfo:Formal class 'renderInfo' [package "graph"] with 4 slots
  .. .. .. ..@ nodes: list()
  .. .. .. ..@ edges: list()
  .. .. .. ..@ graph: list()
  .. .. .. ..@ pars : list()
  .. ..@ graphData :List of 1
  .. .. ..$ edgemode: chr "directed"
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
data(gmI)
## maybe str(gmI) ; plot(gmI) ...
}
\keyword{datasets}
