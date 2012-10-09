#' Tests the calculation of BIC and MLE
#' 
#' @author Alain Hauser
#' $Id: $

cat("Testing the calculation of the BIC score... ")

library(pcalg, lib.loc = "~/R/2.15.1/library")
library(Matrix)

load("test_bicscore.rda") # in directory tests/ i.e., typically *not* installed

score <- new("gauss.l0pen.int.score", 
    targets = gauss.targets, 
    target.index = gauss.target.index, 
    data = gauss.data)

# print(score$pp.dat)

if (any(score$pp.dat$data.count != 1000))
  stop("The number of non-interventions are not calculated correctly.")

if (any(score$pp.dat$scatter.index != 1:5))
  stop("The indices of the scatter matrices are not calculated correctly.")

for (i in 1:5)
  if (norm(score$pp.dat$scatter[[score$pp.dat$scatter.index[i]]][1:5, 1:5] - gauss.scatter[[i]]) > 1e-6)
    stop("The scatter matrices are not calculated correctly.")

for (i in 1:5)
  if (gauss.loc.score[i] != score$local.score(i, gauss.parents[[i]]))
    stop("The local score is not calculated correctly.")

cat("Done.")
