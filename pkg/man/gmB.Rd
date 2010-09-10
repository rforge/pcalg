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
  .. ..@ edgeData  :Formal class 'attrData' [package "graph"] with 2 slots
  .. ..@ nodeData  :Formal class 'attrData' [package "graph"] with 2 slots
  .. ..@ renderInfo:Formal class 'renderInfo' [package "graph"] with 4 slots
  .. ..@ graphData :List of 1
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
% \source{
% %%  ~~ reference to a publication or URL from which the data were obtained ~~
% }
% \references{
% %%  ~~ possibly secondary sources and usages ~~
% }
\examples{
data(gmB)
## maybe str(gmB) ; plot(gmB) ...
}
\keyword{datasets}
