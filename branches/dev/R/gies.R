## GIES algorithm
##
## Author: Alain Hauser <alain.hauser@bfh.ch>
## $Id$
###############################################################################

##################################################
## Auxiliary functions for simulations
##################################################

#' Randomly generates a Gaussian causal model >>>  ../man/r.gauss.pardag.Rd
#'
#' @param p number of vertices
#' @param prob probability of inserting an edge between two given
#'                    vertices
#' @param top.sort indicates whether the produced DAG should be
#'                    topologically sorted
#' @param normalize indicates whether weights and error variances
#'                    should be normalized s.t. the diagonal of the
#'                    corresponding covariance matrix is 1. Note that
#'                    weights and error variances can then lie outside
#'                    the boundaries specified below!
#' @param lbe lower bound of edge weights. Default: 0.1
#' @param ube upper bound of edge weights. Default: 1
#' @param neg.coef indicates whether also negative edge weights should
#'                    be sampled
#' @param labels
#' @param lbv lower bound of vertex variance. Default: 0.5
#' @param ubv upper bound of vertex variance. Default: 1
#' @return  an instance of gauss.pardag
r.gauss.pardag <- function(p,
    prob,
    top.sort = FALSE,
    normalize = FALSE,
    lbe = 0.1,
    ube = 1,
    neg.coef = TRUE,
    labels = as.character(1:p),
    lbv = 0.5,
    ubv = 1)
{
  ## Error checking
  stopifnot(is.numeric(p), length(p) == 1, p >= 2,
      is.numeric(prob), length(prob) == 1, 0 <= prob, prob <= 1,
      is.numeric(lbe), is.numeric(ube), lbe <= ube,
      is.logical(neg.coef),
      is.numeric(lbv), is.numeric(ubv), lbv <= ubv,
      is.character(labels), length(labels) == p)

  ## Create list of nodes, edges and parameters
  edL <- as.list(labels)
  names(edL) <- labels

  ## Create list of parameters; first entry: error variances
  pars <- as.list(runif(p, min = lbv, max = ubv))
  names(pars) <- labels

  ## Create topological ordering
  top.ord <- if (top.sort) 1:p else sample.int(p)

  ## Sample edges and corresponding coefficients, respecting the generated
  ## topological ordering
  for (i in 2:p) {
    ii <- top.ord[i]
    parentCount <- rbinom(1, i - 1, prob)
    edL[[ii]] <- top.ord[sample.int(i - 1, size = parentCount)]
    weights <- runif(parentCount, min = lbe, max = ube)
    if (neg.coef)
      weights <- weights * sample(c(-1, 1), parentCount, replace = TRUE)
    pars[[ii]] <- c(pars[[ii]], 0, weights)
  }
  edL[[top.ord[1]]] <- integer(0)
  pars[[top.ord[1]]] <- c(pars[[top.ord[1]]], 0)

  ## Create new instance of gauss.pardag
  result <- new("GaussParDAG", nodes = labels, in.edges = edL, params = pars)

  ## Normalize if requested
  if (normalize) {
    H <- diag(result$cov.mat())
    result$set.err.var(result$err.var() / H)
    H <- sqrt(H)
    for (i in 1:p)
      if (length(edL[[i]]) > 0)
        result$.params[[i]][-c(1, 2)] <- pars[[i]][-c(1, 2)] * H[edL[[i]]] / H[i]
  }

  ## Validate and return object
  validObject(result)
  result
}

#' Simulates independent observational or interventional data for a
#' specified interventions from a Gaussian causal model
#'
#' @param   n         number of data samples
#' @param   object    an instance of gauss.pardag
#' @param   target    intervention target
#' @param   target.value    value of intervention targets
rmvnorm.ivent <- function(n, object, target = integer(0), target.value = numeric(0))
{

  p <- object$node.count()
  ## Error checking
  stopifnot(length(target) == 0 || (1 <= min(target) && max(target) <= p))
  stopifnot((is.vector(target.value) && length(target.value) == length(target)) ||
            (is.matrix(target.value) && dim(target.value) == c(n, length(target))))

  ## Simulate error terms
  sigma <- sqrt(object$err.var())
  mu <- object$intercept()
  Y <- matrix(rnorm(n*p, mu, sigma), nrow = p, ncol = n)

  ## Insert intervention values
  Y[target, ] <- target.value

  ## Calculate matrix of structural equation system
  A <- - t(object$weight.mat(target))
  diag(A) <- 1.

  ## Solve linear structural equations
  t(solve(A, Y))
}


##################################################
## Structure learning algorithms
##################################################

##' Wrapper function for all causal inference algorithms.  It's not recommended
##' to use it directly; adapted wrapper functions for the single algorithms are
##' provided
#'
##' @param algorithm 	name of the causal inference algorithm to be used
##' @param score 	scoring object to be used
##' @param labels 	node labels
##' @param targets 	unique list of targets. Normally determined from the scoring object
##' @param ... 		additional parameters passed to the algorithm chosen
caus.inf <- function(algorithm = c("GIES", "GDS", "SiMy"),
                     score, labels = score$getNodes(), targets = score$getTargets(), ...)
{
  algorithm <- match.arg(algorithm)
  essgraph <- new("EssGraph", nodes = labels, targets = targets, score = score)
  if (essgraph$caus.inf(algorithm, ...)) {
    if (algorithm == "GIES") {
      ## GIES yields an essential graph; calculate a representative thereof
      list(essgraph = essgraph, repr = essgraph$repr())
    } else {
      ## GDS and SiMy yield a DAG; calculate the corresponding essential graph,
      ## although calculations may come from a model class where Markov equivalence
      ## does not hold!
      list(essgraph = dag2essgraph(essgraph$repr(), targets = targets),
           repr = essgraph$repr())
    }
  } else stop("invalid 'algorithm' or \"EssGraph\" object")
}

##' Greedy Interventional Equivalence Search - GIES --> ../man/gies.Rd
##'
##' @param score	scoring object to be used
##' @param labels	node labels
##' @param targets	unique list of targets. Normally determined from the scoring object
##' @param fixedGaps	logical matrix indicating forbidden edges
##' @param adaptive sets the behaviour for adaptiveness in the forward phase (cf. "ARGES")
##' @param turning	indicates whether the turning step should be indicated.
##' @param maxDegree	maximum vertex degree allowed
##' @param verbose	indicates whether debug output should be printed
##' @param ...		additional parameters (currently none)
gies <- function(
    score, 
    labels = score$getNodes(), 
    targets = score$getTargets(),
    fixedGaps = NULL, 
    adaptive = c("none", "vstructures", "triples"), 
    turning = TRUE, 
    maxDegree = integer(0),
    verbose = FALSE, 
    ...)
{
  caus.inf("GIES", score = score, labels = labels, targets = targets,
           fixedGaps = fixedGaps, adaptive = adaptive, turning = turning,
           maxDegree = maxDegree, verbose = verbose, ...)
}

##' Greedy Equivalence Search - GES --> ../man/ges.Rd
##'
##' @param score 	scoring object to be used
##' @param labels 	node labels
##' @param fixedGaps 	logical matrix indicating forbidden edges
##' @param adaptive sets the behaviour for adaptiveness in the forward phase (cf. "ARGES")
##' @param turning 	indicates whether the turning step should be indicated.
##' 		Setting this parameter to FALSE gives Chickering's original version
##' @param maxDegree 	maximum vertex degree allowed
##' @param verbose 	indicates whether debug output should be printed
##' @param ... 		additional parameters (currently none)
##' @param targets 	unique list of targets. Normally determined from the scoring object
ges <- function(
    score, 
    labels = score$getNodes(),
    fixedGaps = NULL, 
    adaptive = c("none", "vstructures", "triples"), 
    turning = TRUE, 
    maxDegree = integer(0),
    verbose = FALSE, 
    ...)
{
  caus.inf("GIES", score = score, labels = labels, targets = list(integer(0)),
           fixedGaps = fixedGaps, adaptive = adaptive, turning = turning,
           maxDegree = maxDegree, verbose = verbose, ...)
}

##' Greedy DAG Search - GDS : greedy search in the DAG space --> ../man/gds.Rd
##'
##' @param score 	scoring object to be used
##' @param labels 	node labels
##' @param targets
##' @param fixedGaps 	logical matrix indicating forbidden edges
##' @param turning 	indicates whether the turning step should be indicated.
##' 		Setting this parameter to FALSE gives Chickering's original version
##' @param maxDegree 	maximum vertex degree allowed
##' @param verbose 	indicates whether debug output should be printed
##' @param ... 		additional parameters (currently none)
gds <- function(
    score, 
    labels = score$getNodes(), 
    targets = score$getTargets(),
    fixedGaps = NULL, 
    turning = TRUE, 
    maxDegree = integer(0), 
    verbose = FALSE, 
    ...)
{
  caus.inf("GDS", score = score, labels = labels, targets = targets,
           fixedGaps = fixedGaps, turning = turning, maxDegree = maxDegree, verbose = verbose, ...)
}

##' Dynamic programming approach of Silander and Myllimäki - SiMy --> ../man/simy.Rd
##'
##' @param score 	scoring object to be used
##' @param labels 	node labels
##' @param targets
##' @param verbose 	indicates whether debug output should be printed
##' @param ... 		additional parameters (currently none)
simy <- function(score, labels = score$getNodes(), targets = score$getTargets(),
                 verbose = FALSE, ...)
{
  caus.inf("SiMy", score = score, labels = labels, targets = targets, verbose = verbose, ...)
}


#' Converts a DAG to an (observational or interventional) essential graph
dag2essgraph <- function(dag, targets = list(integer(0))) {
  edgeListDAG <- inEdgeList(dag)
  edgeListEssGraph <- .Call("dagToEssentialGraph", edgeListDAG, targets)
  if (is.matrix(dag)) {
    p <- nrow(dag)
    result <- sapply(1:p, function(i) 1:p %in% edgeListEssGraph[[i]])
    rownames(result) <- rownames(dag)
    colnames(result) <- colnames(dag)
    result
  } else if (inherits(dag, "graphNEL")) {
    nodeNames <- nodes(dag)
    names(edgeListEssGraph) <- nodeNames
    result <- new("graphNEL",
        nodes = nodeNames,
        edgeL = lapply(edgeListEssGraph, function(v) nodeNames[v]),
        edgemode = "directed")
    reverseEdgeDirections(result)
  } else {
    new("EssGraph",
        nodes = dag$.nodes,
        in.edges = edgeListEssGraph,
        targets = targets)
  }
}

