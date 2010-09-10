\name{gmG}
\alias{gmG}
\docType{data}
\title{Graphical Model 8-Dimensional Gaussian Example Data}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(gmG)}
\format{
  The format is:
List of 2
 $ x: num [1:5000, 1:8] 0.181 0.354 0.504 0.701 -0.705 ...
 $ g:Formal class 'graphNEL' [package "graph"] with 6 slots
  .. ..@ nodes     : chr [1:8] "1" "2" "3" "4" ...
  .. ..@ edgeL     :List of 8
  .. .. ..$ 1:List of 1
  .. .. .. ..$ edges: int [1:3] 8 6 2
  .. .. ..$ 2:List of 1
  .. .. .. ..$ edges: int [1:2] 5 3
  .. .. ..$ 3:List of 1
  .. .. .. ..$ edges: int(0)
  .. .. ..$ 4:List of 1
  .. .. .. ..$ edges: int(0)
  .. .. ..$ 5:List of 1
  .. .. .. ..$ edges: int [1:2] 8 6
  .. .. ..$ 6:List of 1
  .. .. .. ..$ edges: int 7
  .. .. ..$ 7:List of 1
  .. .. .. ..$ edges: int(0)
  .. .. ..$ 8:List of 1
  .. .. .. ..$ edges: int(0)
  .. ..@ edgeData  :Formal class 'attrData' [package "graph"] with 2 slots
  .. .. .. ..@ data    :List of 8
  .. .. .. .. ..$ 1|8:List of 1
  .. .. .. .. .. ..$ weight: num 0.276
  .. .. .. .. ..$ 1|6:List of 1
  .. .. .. .. .. ..$ weight: num 0.515
  .. .. .. .. ..$ 1|2:List of 1
  .. .. .. .. .. ..$ weight: num 0.283
  .. .. .. .. ..$ 2|5:List of 1
  .. .. .. .. .. ..$ weight: num 0.187
  .. .. .. .. ..$ 2|3:List of 1
  .. .. .. .. .. ..$ weight: num 0.733
  .. .. .. .. ..$ 5|8:List of 1
  .. .. .. .. .. ..$ weight: num 0.979
  .. .. .. .. ..$ 5|6:List of 1
  .. .. .. .. .. ..$ weight: num 0.186
  .. .. .. .. ..$ 6|7:List of 1
  .. .. .. .. .. ..$ weight: num 0.372
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
data(gmG)
## maybe str(gmG) ; plot(gmG) ...
}
\keyword{datasets}
