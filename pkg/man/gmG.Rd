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
data(gmG)
## maybe str(gmG) ; plot(gmG) ...
}
\keyword{datasets}
