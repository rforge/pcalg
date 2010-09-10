\name{gmD}
\alias{gmD}
\docType{data}
\title{Graphical Model Discrete 5-Dim Example Data}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(gmD)}
\format{
  The format is:
List of 2
 $ x:'data.frame':	10000 obs. of  5 variables:
  ..$ X1: int [1:10000] 2 2 1 1 1 2 2 0 2 0 ...
  ..$ X2: int [1:10000] 0 1 0 0 0 1 1 0 0 0 ...
  ..$ X3: int [1:10000] 0 1 1 2 0 1 0 1 2 0 ...
  ..$ X4: int [1:10000] 2 2 3 2 2 0 2 2 2 0 ...
  ..$ X5: int [1:10000] 1 1 0 1 1 1 1 1 1 0 ...
 $ g:Formal class 'graphNEL' [package "graph"] with 6 slots
  .. ..@ nodes     : chr [1:5] "1" "2" "3" "4" ...
  .. ..@ edgeL     :List of 5
  .. .. ..$ 1:List of 1
  .. .. .. ..$ edges: int 3
  .. .. ..$ 2:List of 1
  .. .. .. ..$ edges: int [1:2] 3 4
  .. .. ..$ 3:List of 1
  .. .. .. ..$ edges: int(0)
  .. .. ..$ 4:List of 1
  .. .. .. ..$ edges: int 5
  .. .. ..$ 5:List of 1
  .. .. .. ..$ edges: int(0)
  .. ..@ edgeData  :Formal class 'attrData' [package "graph"] with 2 slots
  .. .. .. ..@ data    :List of 4
  .. .. .. .. ..$ 1|3:List of 1
  .. .. .. .. .. ..$ weight: num 1
  .. .. .. .. ..$ 2|3:List of 1
  .. .. .. .. .. ..$ weight: num 1
  .. .. .. .. ..$ 2|4:List of 1
  .. .. .. .. .. ..$ weight: num 1
  .. .. .. .. ..$ 4|5:List of 1
  .. .. .. .. .. ..$ weight: num 1
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
data(gmD)
## maybe str(gmD) ; plot(gmD) ...
}
\keyword{datasets}
