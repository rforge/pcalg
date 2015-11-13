#' Tests the calculation of BIC and MLE
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

for (m in 1:nrow(settings)) {
  cat(sprintf("Setting: storage format = %s, C++ library = %s\n", 
          settings$format[m], settings$cpp[m]))
  
  score <- new("GaussL0penIntScore", 
      targets = gauss.targets, 
      target.index = gauss.target.index, 
      data = gauss.data,
      format = settings$format[m],
      use.cpp = settings$cpp[m],
      intercept = FALSE)
  
  # print(score$pp.dat)

  if (any(score$pp.dat$data.count != 1000))
    stop("The number of non-interventions are not calculated correctly.")

  if (settings$format[m] == "scatter") {
    if (any(score$pp.dat$scatter.index != 1:5))
      stop("The indices of the scatter matrices are not calculated correctly.")
  
    for (i in 1:5)
      if (!isTRUE(all.equal(score$pp.dat$scatter[[score$pp.dat$scatter.index[i]]][1:5, 1:5], 
          gauss.scatter[[i]], 
          tolerance = tol)))
        stop("The scatter matrices are not calculated correctly.")
  } # IF "scatter"
  
  for (i in 1:5) {
    if (!isTRUE(all.equal(gauss.loc.score[[i]], 
        score$local.score(i, gauss.parents[[i]]), 
        tolerance = tol)))
      stop("The local score is not calculated correctly.")
  }
  
  # print(lapply(1:5, function(i) score$local.fit(i, gauss.parents[[i]])))
  
  for (i in 1:5) {
    local.mle <- score$local.fit(i, gauss.parents[[i]])
    if (length(local.mle) != length(gauss.mle[[i]]) ||
        !isTRUE(all.equal(gauss.mle[[i]], local.mle, 
        tolerance = tol)))
      stop("The local MLE is not calculated correctly.")
  }
}

cat("Done.\n")
