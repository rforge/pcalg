\name{legal.path}
\alias{legal.path}
\title{
Check if a 3-node-path is legal.
}
\description{
  Check if the path a-b-c is legal. I.e., either b is a collider or
  a-b-c is a triangle.
}
\usage{
legal.path(a, b, c, amat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a}{Node a.
}
  \item{b}{Node b.
}
  \item{c}{Node c.
}
  \item{amat}{Adjecency matrix (coding 0,1,2 for no edge, circle, arrowhead)
}
}
\value{\code{TRUE} if path is legal, otherwise \code{FALSE}}
}
\author{
  Markus Kalisch (\email{kalisch@stat.math.ethz.ch})
}
\examples{
amat <- matrix( c(0,1,1,0,0, 2,0,1,0,0, 2,2,0,2,1, 0,0,1,0,0,
0,0,2,0,0), 5,5)

legal.path(1,3,5, amat)
legal.path(1,2,3, amat)
legal.path(2,3,4, amat)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}
