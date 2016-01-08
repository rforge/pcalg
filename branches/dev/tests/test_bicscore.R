#' Tests the calculation of BIC and MLE, as well as basic functions
#' the corresponding score class
#' 
#' @author Alain Hauser
#' $Id$

cat("Testing the calculation of the BIC score...\n")

library(pcalg)

load("test_bicscore.rda") # in directory tests/ i.e., typically *not* installed

## Tolerance for numerical comparison
tol <- sqrt(.Machine$double.eps)

## Define all test settings
settings <- expand.grid(
  format = c("scatter", "raw"),
  cpp = c(FALSE, TRUE),
  stringsAsFactors = FALSE)
nreps <- 5

for (m in 1:nrow(settings)) {
  cat(sprintf("Setting: storage format = %s, C++ library = %s\n", 
          settings$format[m], settings$cpp[m]))
  
  for (i in 1:nreps) {
    perm <- 1:nrow(gauss.data)
    
    ## Randomly permute data
    if (i > 1) {
      set.seed(i)
      perm <- sample(perm)
    }
    
    ## Try to create the score object with non-valid data,
    ## check if error is thrown
    
    
    ## Create the score object with valid data
    score <- new("GaussL0penIntScore", 
        targets = gauss.targets, 
        target.index = gauss.target.index[perm], 
        data = gauss.data[perm, ],
        format = settings$format[m],
        use.cpp = settings$cpp[m],
        intercept = FALSE)
    
    # print(score$pp.dat)
  
    if (any(score$pp.dat$data.count != 1000))
      stop("The number of non-interventions are not calculated correctly.")
  
    if (settings$format[m] == "scatter") {
      if (any(score$pp.dat$scatter.index != 1:5))
        stop("The indices of the scatter matrices are not calculated correctly.")
    
      for (j in 1:5)
        if (!isTRUE(all.equal(score$pp.dat$scatter[[score$pp.dat$scatter.index[j]]][1:5, 1:5], 
            gauss.scatter[[j]], 
            tolerance = tol)))
          stop("The scatter matrices are not calculated correctly.")
    } # IF "scatter"
    
    for (j in 1:5) {
      if (!isTRUE(all.equal(gauss.loc.score[[j]], 
          score$local.score(j, gauss.parents[[j]]), 
          tolerance = tol)))
        stop("The local score is not calculated correctly.")
    }
    
    # print(lapply(1:5, function(i) score$local.fit(i, gauss.parents[[i]])))
    
    for (j in 1:5) {
      local.mle <- score$local.fit(j, gauss.parents[[j]])
      if (length(local.mle) != length(gauss.mle[[j]]) ||
          !isTRUE(all.equal(gauss.mle[[j]], local.mle, 
          tolerance = tol)))
        stop("The local MLE is not calculated correctly.")
    }
  }
}

cat("Done.\n")
