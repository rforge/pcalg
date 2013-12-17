\name{trueCov}
\alias{trueCov}
\title{
Covariance matrix of a DAG.
}
\description{
Compute the (true) covariance matrix of a generated DAG.
}
\usage{
trueCov(dag)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dag}{Graph object containing the DAG.}
}
\value{
Covariance matrix.
}
\author{
Markus Kalisch
}
\note{
This function can \emph{not} be used to estimate the covariance matrix
from an estimated DAG or corresponding data.
}
\seealso{
  \code{\link{randomDAG}} for generating a random DAG
}
\examples{
set.seed(123)
g <- randomDAG(n = 5, prob = 0.3) ## generate random DAG
if(require(Rgraphviz)) {
plot(g)
}

## Compute true covariance matrix
trueCov(g)

## For comparison:
## Estimate true covariance matrix after generating data from g
d <- rmvDAG(10000, g)
cov(d)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}
