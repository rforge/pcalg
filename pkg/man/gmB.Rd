\name{gmB}
\alias{gmB}
\docType{data}
\title{Graphical Model 5-Dim Binary Example Data}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(gmB)}
\format{
  The format is:
List of 2
 $ x: int [1:5000, 1:5] 0 1 1 0 0 1 1 0 1 1 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : NULL
  .. ..$ : chr [1:5] "V1" "V2" "V3" "V4" ...
 $ g:Formal class 'graphNEL' [package "graph"] with 6 slots
  .. ..@ nodes     : chr [1:5] "1" "2" "3" "4" ...
  .. ..@ edgeL     :List of 5
  .. .. ..$ 1:List of 1
  .. .. .. ..$ edges: int(0)
  .. .. ..$ 2:List of 1
  .. .. .. ..$ edges: int 1
  .. .. ..$ 3:List of 1
  .. .. .. ..$ edges: int 5
  .. .. ..$ 4:List of 1
  .. .. .. ..$ edges: int [1:2] 2 3
  .. .. ..$ 5:List of 1
  .. .. .. ..$ edges: int 2
  .. ..@ edgeData  :Formal class 'attrData' [package "graph"] with 2 slots
  .. .. .. ..@ data    :List of 5
  .. .. .. .. ..$ 2|1:List of 1
  .. .. .. .. .. ..$ weight: num 1
  .. .. .. .. ..$ 3|5:List of 1
  .. .. .. .. .. ..$ weight: num 1
  .. .. .. .. ..$ 4|2:List of 1
  .. .. .. .. .. ..$ weight: num 1
  .. .. .. .. ..$ 4|3:List of 1
  .. .. .. .. .. ..$ weight: num 1
  .. .. .. .. ..$ 5|2:List of 1
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
data(gmB)
## maybe str(gmB) ; plot(gmB) ...
}
\keyword{datasets}
