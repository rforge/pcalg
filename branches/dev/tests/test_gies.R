#' Tests the causal inference algorithms for interventional data:
#' GIES, GES, DP
#' 
#' @author Alain Hauser
#' $Id$

cat("Testing the causal inference algorithms for interventional data: ")

library(pcalg)

load("test_bicscore.rda") # in directory tests/ i.e., typically *not* installed
p <- ncol(gauss.data)

## Tolerance for numerical comparison
tol <- sqrt(.Machine$double.eps)

fcns <- c(gies, gds)
fcn.names <- c("GIES", "GDS")

for (fi in 1:length(fcns)) {
  cat(paste(fcn.names[fi], "... ", sep = ""))
  
  for (cpp in c(FALSE, TRUE)) {
    score <- new("GaussL0penIntScore", 
      targets = gauss.targets, 
      target.index = gauss.target.index, 
      data = gauss.data,
      use.cpp = cpp)
  
    # gies(p, gauss.targets, score, DEBUG.LEVEL = 3)
    est.graph <- fcns[[fi]](p, gauss.targets, score)
    
    for (i in 1:p)
      if (!isTRUE(all.equal(est.graph$essgraph$.in.edges[[i]], gauss.parents[[i]], tolerance = tol)))
        stop("Parents are not estimated correctly.")
  }
}

cat("Done.")
