#' Tests the causal inference algorithms for interventional data:
#' GIES, GES, DP
#' 
#' @author Alain Hauser
#' $Id$

cat("Testing the causal inference algorithms for interventional data:\n")

library(pcalg)

load("test_bicscore.rda") # in directory tests/ i.e., typically *not* installed
str(gauss.data)
p <- ncol(gauss.data)

(doExtras <- pcalg:::doExtras())
DBG <- if(doExtras) 2 else 0 # no debugging by default
## Tolerance for numerical comparison
tol <- sqrt(.Machine$double.eps) # = default for all.equal()

fcns <- c(GIES = gies, GDS = gds)

for (nf in names(fcns)) {
  cat(if(doExtras)"\n\n", nf, if(doExtras)":\n" else ": ... ",
      if(doExtras) paste0(paste(rep("=", nchar(nf)), collapse=""), "\n"),
      sep = "")
  for (cpp in c(FALSE, TRUE)) {
    score <- new("GaussL0penIntScore", 
                 targets = gauss.targets, 
                 target.index = gauss.target.index, 
                 data = gauss.data,
                 use.cpp = cpp)
    est.graph <- fcns[[nf]](p, gauss.targets, score, DEBUG.LEVEL = DBG)
    for (i in 1:p) {
      if(doExtras) cat("  use.cpp = ", cpp,"; i = ", i, "\n", sep="")
      if (!isTRUE(all.equal(est.graph$essgraph$.in.edges[[i]],
                            gauss.parents[[i]], tolerance = tol)))
        stop("Parents are not estimated correctly.")
    }
    print(proc.time())
  }
  cat("[Ok]\n")
}
cat(if(doExtras) "\n", "Done.\n")
