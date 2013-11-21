library(pcalg)

## NB: add tests in addition to the simple onesfrom Maathuis and Colombo (2013)
## in ../man/backdoor.Rd

`%w/o%` <- function(x, y) x[!x %in% y] #--  x without y
## slightly faster:
`%w/o%` <- function(x, y) x[!match(x, y, nomatch = 0L)]

set.seed(47)
p <- 17
myDAG <- randomDAG(p, prob = 1/4) ## true DAG

## Extract the adjacency matrix of the true DAG
true.amat <- (amat <- as(myDAG, "matrix")) != 0 # TRUE/FALSE <==> 1/0
print.table(1*true.amat, zero.=".") # "visualization"

nodes <- 1:p; names(nodes) <- nodes
cat("Time for many backdoor() s : ", system.time(
LL <- lapply(nodes, function(i)
	     lapply(nodes %w/o% i,
		    backdoor,
		    amat = true.amat, x = i, type="dag"))
), "\n")

for(i in nodes[1:3]) ## Nodes 1,2,3 are all "root" nodes:
    stopifnot(vapply(LL[[i]], identical, NA, y=integer(0)))

str(LL[-(1:3)]) ## Martin: interesting.. is "this" known?



