##################################################
## Classes
##################################################
setClass("gAlgo",
         representation(call = "call",
                        n	   = "integer",
                        max.ord = "integer",
                        n.edgetests= "numeric",
                        sepset= "list",
                        pMax= "matrix"), "VIRTUAL")


setClass("fciAlgo",
         representation(amat = "matrix", allPdsep = "list",
                        n.edgetestsPDSEP = "numeric", max.ordPDSEP = "integer"),
         contains = "gAlgo")

setClass("pcAlgo",
         representation(graph = "graph", zMin = "matrix"), ## zMin for compatibility
         contains = "gAlgo")

##################################################
## Reference classes used by GIES
##################################################

#' Virtual base class for all parametric causal models.
#' The meaning of the "params" depends on the model used.
setRefClass("pardag",
    fields = list(
        .nodes = "vector",
        .in.edges = "list",
        .params = "list"),
    
    validity = function(object) {
      if (anyDuplicated(object$.nodes))
        return("The node names must be unique")
      if (any(names(object$.in.edges) != object$.nodes))
        return("The elements of 'in.edges' must be named after the nodes.")
      if (!all(sapply(object$.in.edges, is.numeric)))
        return("The vectors in 'in.edges' must contain numbers.")
      
      edgeRange <- range(unlist(object$.in.edges))
      if (object$edge.count() > 0 &&
          (edgeRange[1] < 1 || edgeRange[2] > object$node.count()))
        return("Invalid range of edge sources.")
      
      return(TRUE)
    },
    
    methods = list(
        #' Constructor
        initialize = function(nodes, in.edges = NULL, params = list()) {
          .nodes <<- nodes
          
          if (is.null(in.edges))
            .in.edges <<- replicate(length(nodes), integer(0), simplify = FALSE)
          else
            .in.edges <<- lapply(1:length(in.edges), function(i) as.integer(in.edges[[i]]))
          names(.in.edges) <<- nodes
          
          .params <<- params
        },
        
        #' Yields the number of nodes
        node.count = function() {
          length(.nodes)
        },
        
        #' Yields the total number of edges in the graph
        edge.count = function() {
          sum(sapply(.in.edges, length))
        }
        ),
        
    "VIRTUAL")

#' Coercion to a graphNEL instance
setAs("pardag", "graphNEL",
    def = function(from) {
      result <- new("graphNEL", 
          nodes = from$.nodes, 
          edgeL = from$.in.edges, 
          edgemode = "directed")
      return(reverseEdgeDirections(result))
    })

#' Coercion to a (logical) matrix
setAs("pardag", "matrix",
    def = function(from) {
      p <- from$node.count()
      sapply(1:p, function(i) 1:p %in% from$.in.edges[[i]])
    })

#' Plot method (needs Rgraphviz to work!!)
setMethod("plot", "pardag", 
    function(x, y, ...) {
      if (!validObject(x))
        stop("The parametric DAG model to be plotted is not valid")
      
      if (missing(y))
        y <- "dot"
      invisible(plot(as(x, "graphNEL"), y, ...))
    })


#' Virtual base class for all scoring classes
setRefClass("int.score",
    fields = list(
        .node.count = "integer",
        .targets = "list",
        .target.index = "vector",
        .data = "matrix",
        decomp = "logical",
        c.fcn = "character",
        .pardag.class = "character"),
    
    validity = function(object) {
      ## Check if targets are valid (i.e., unique)
      targets.tmp <- object$.targets
      for (i in seq(along = targets.tmp)) {
        targets.tmp[[i]] <- sort(unique(targets.tmp[[i]]))
        if (length(targets.tmp[[i]]) != length(object$.targets[[i]]))
          return("Target variables must not be listed multiple times.")
      }
      if (length(unique(targets.tmp)) != length(targets.tmp))
        return("Targets must not be listed multiple times.")
      
      ## Check whether data is available from all intervention targets
      if (unique(object$.target.index) != 1:length(object$.targets))
        return("Data from all intervention targets must be available")
      
      ## Check if dimensions of target.index and data conincide
      if (length(object$.target.index) != nrow(object$.data))
        return("Length of target index vector does not coincide with sample size.")
      
      return(TRUE)
    },
    
    methods = list(
        #' Constructor
        initialize = function(targets, target.index, data, ...) {
          ## Order by ascending target indices (necessary for certain scoring objects)
          if (is.unsorted(target.index)) 
            perm <- order(target.index)
          else
            perm <- 1:length(target.index)
          
          .targets <<- targets
          .target.index <<- target.index[perm]
          .data <<- data[perm, ]
          .node.count <<- ncol(data)
          
          ## Declare scores as not decomposable "by default"
          decomp <<- FALSE
              
          ## No C++ scoring object by default
          c.fcn <<- "none"

          callSuper(...)
        },
                
        #' Calculates the local score of a vertex and its parents
        local.score = function(vertex, parents) {
          stop("local.score is not implemented in this class.")
        },
        
        #' Calculates the global score of a DAG
        global.score = function(dag) {
          sum(sapply(1:.node.count,
                  function(i) local.score(i, dag$.in.edges[[i]])))
        },
        
        #' Calculates the global score of a DAG which is only specified
        #' by its list of in-edges
        global.score.int = function(edges) {
          sum(sapply(1:.node.count,
                  function(i) local.score(i, edges[[i]])))
        },
        
        #' Calculates the local MLE for a vertex and its parents
        local.mle = function(vertex, parents) {
          stop("local.mle is not implemented in this class.")
        },
        
        #' Calculates the global MLE
        global.mle = function(dag) {
          lapply(1:.node.count,
              function(i) local.mle(i, dag$.in.edges[[i]]))
        }
        ),
    "VIRTUAL")

#' BIC score for Gaussian causal models
setRefClass("gauss.bic.int.score",
    contains = "int.score",
    
    fields = list(
        .intercept = "logical",
        .scatter = "list",
        .scatter.index = "vector",
        .data.count = "vector",
        .total.data.count = "integer"),
    
    validity = function(object) {
      if (unique(object$.scatter.index) != 1:length(object$.scatter))
        return("The index list of distinct scatter matrices has an invalid range.")
      p <- ncol(object$.data)
      if (any(sapply(1:length(object$.scatter),
              function(i) dim(object$.scatter[[i]]) != c(p, p))))
        return("The scatter matrices have invalid dimensions.")
      
      return(TRUE)
    },
    
    methods = list(
        #' Constructor
        initialize = function(targets, target.index, data, intercept = FALSE, ...) {
          ## Store supplied data in sorted form
          callSuper(targets, target.index, data)
          
          ## BIC is decomposable
          decomp <<- TRUE
          
          ## Underlying causal model class: Gaussian
          .pardag.class <<- "gauss.pardag"
          
          ## If an intercept is allowed, add column of ones to data matrix
          .intercept <<- intercept
          if (intercept)
            data <- cbind(data, 1)
          
          ## Number of variables
          p <- ncol(data)
          .total.data.count <<- as.integer(nrow(data))
          
          ## Create scatter matrices for different targets
          ti.lb <- c(sapply(1:length(.targets), function(i) match(i, .target.index)), 
              length(.target.index) + 1)
          scatter.mat <- lapply(1:length(.targets), 
              function(i) crossprod(data[ti.lb[i]:(ti.lb[i + 1] - 1), ]))
          
          ## Find all interventions in which the different variables
          ## are _not_ intervened
          non.ivent <- matrix(FALSE, ncol = p, nrow = length(.targets))
          .scatter.index <<- integer(p)
          .data.count <<- integer(p)
          max.si <- 0
          for (i in 1:p) {
            ## Generate indices of (distinct) scatter matrices
            non.ivent[ , i] <- sapply(1:length(.targets), 
                function(j) !(i %in% .targets[[j]]))
            .scatter.index[i] <<- max.si + 1
            j <- 1
            while (j < i) {
              if (all(non.ivent[, i] == non.ivent[, j])) {
                .scatter.index[i] <<- .scatter.index[j]
                j <- i
              }
              j <- j + 1
            }
            if (.scatter.index[i] == max.si + 1)
              max.si <- max.si + 1
            
            ## Count data samples from "non-interventions" at i
            .data.count[i] <<- sum(ti.lb[which(non.ivent[, i]) + 1] - ti.lb[which(non.ivent[, i])])
          }
                    
          ## Calculate the distinct scatter matrices for the
          ## "non-interventions"
          .scatter <<- lapply(1:max.si, 
              function(i) Reduce("+", scatter.mat[non.ivent[, match(i, .scatter.index)]]))
        },
        
        #' Calculates the local score of a vertex and its parents
        local.score = function(vertex, parents) {
          ## If an intercept is allowed, add a fake parent node
          parents <- sort(parents)
          if (.intercept)
            parents <- c(.node.count + 1, parents)
          
          sigma2 <- .scatter[[.scatter.index[vertex]]][vertex, vertex]
          if (length(parents) != 0) {
            b <- .scatter[[.scatter.index[vertex]]][vertex, parents]
            sigma2 <- sigma2 - as.numeric(b %*% solve(.scatter[[.scatter.index[vertex]]][parents, parents], b))
          }
          
          -0.5*(.data.count[vertex]*(1 + log(sigma2/.data.count[vertex])) + log(.total.data.count)*(1 + length(parents)))
        },
        
        #' Calculates the local MLE for a vertex and its parents
        local.mle = function(vertex, parents) {
          ## If an intercept is allowed, add a fake parent node
          parents <- sort(parents)
          if (.intercept)
            parents <- c(.node.count + 1, parents)
          
          sigma2 <- .scatter[[.scatter.index[vertex]]][vertex, vertex]
          if (length(parents) != 0) {
            beta <- solve(.scatter[[.scatter.index[vertex]]][parents, parents], 
                .scatter[[.scatter.index[vertex]]][vertex, parents])
            sigma2 <- sigma2 - .scatter[[.scatter.index[vertex]]][vertex, parents] %*% beta
          }
          else
            beta <- numeric(0)
          
          if (.intercept)
            return(c(sigma2/.data.count[vertex], beta))
          else
            return(c(sigma2/.data.count[vertex], 0, beta))
        }
        )
    )

#' Interventional essential graph
setRefClass("ess.graph",
    fields = list(
        .nodes = "vector",
        .in.edges = "list",
        targets = "list",
        score = "int.score"
    ),
    
    validity = function(object) {
      if (any(names(object$.in.edges) != object$.nodes))
        return("The elements of 'in.edges' must be named after the nodes.")
      if (!all(sapply(object$.in.edges, is.numeric)))
        return("The vectors in 'in.edges' must contain numbers.")
      
      edgeRange <- range(unlist(object$.in.edges))
      if (object$edge.count() > 0 &&
          (edgeRange[1] < 1 || edgeRange[2] > object$node.count()))
        return("Invalid range of edge sources.")
      
      return(TRUE)
    },
    
    methods = list(
        #' Constructor
        initialize = function(nodes, in.edges, ...) {
          .nodes <<- nodes
          if (missing(in.edges))
            .in.edges <<- replicate(length(nodes), integer(0))
          else
            .in.edges <<- in.edges
          names(.in.edges) <<- nodes
          
          callSuper(...)
        },

        #' Yields the number of nodes
        node.count = function() {
          length(.nodes)
        },
        
        #' Yields the total number of edges in the graph
        edge.count = function() {
          sum(sapply(.in.edges, length))
        },
        
        #' Creates a list of options for the C++ function "causalInference";
        #' internal function
        causal.inf.options = function(caching = TRUE, 
            turning = TRUE, 
            maxdegree = integer(0),
            maxsteps = 0,
            childrenonly = integer(0),
            DEBUG.LEVEL = 0) {
          list(caching = caching,
              turning = turning,
              maxdegree = maxdegree,
              maxsteps = maxsteps,
              childrenonly = childrenonly,
              DEBUG.LEVEL = DEBUG.LEVEL)
        },
        
        #' Performs one greedy step
        greedy.step = function(direction = c("forward", "backward", "turning")) {
          ## Cast direction
          direction <- match.arg(direction)
          alg.name <- switch(direction,
              forward = "GIES-F",
              backward = "GIES-B",
              turning = "GIES-T")
                    
          score.fcn <- ifelse(score$decomp, 
              function(vertex, parents) score$local.score(vertex, parents),
              function(edges) score$global.score.int(edges))
          substr(direction, 1, 1) <- toupper(substr(direction, 1, 1))
          
          new.graph <- .Call("causalInference", 
              .in.edges,
              targets,
              alg.name,
              score.fcn, 
              causal.inf.options(caching = FALSE, maxsteps = 1),
              PACKAGE = "pcalg")
          if (new.graph$steps > 0) {
            .in.edges <<- new.graph$.in.edges
            names(.in.edges) <<- .nodes
          }
          
          return(new.graph$steps == 1)
        },
        
        greedy.search = function(direction) {
          score.fcn <- ifelse(score$decomp, 
              function(vertex, parents) score$local.score(vertex, parents),
              function(edges) score$global.score.int(edges))
          substr(direction, 1, 1) <- toupper(substr(direction, 1, 1))
          
          new.graph <- .Call(sprintf("greedy%s", direction), .in.edges, score.fcn)
          if (new.graph$steps > 0)
            .in.edges <<- new.graph$in.edges
          
          return(new.graph$steps)
        },
        
        #' Performs a causal inference from an arbitrary start DAG
        #' with a specified algorithm
        caus.inf = function(algorithm, ...) {
          if (score$decomp)
            score.fcn <- function(vertex, parents) score$local.score(vertex, parents)
          else
            score.fcn <- function(edges) score$global.score.int(edges)
          
          new.graph <- .Call("causalInference", 
              .in.edges,
              targets,
              algorithm,
              score.fcn, 
              causal.inf.options(...),
              PACKAGE = "pcalg")
          
          .in.edges <<- new.graph$in.edges
          names(.in.edges) <<- .nodes
        },
        
        #' Performs GIES from an arbitrary start DAG
        gies = function(...) caus.inf("GIES", ...),
        
        #' Performs GDS from an arbitrary start DAG
        gds = function(...) caus.inf("GDS", ...),
        
        #' DP search of Silander and Myllymäki (ignores the start DAG!)
        silander = function(...) caus.inf("Silander", ...),
        
        #' Yields a representative (estimating parameters via MLE)
        repr = function() {
          in.edges <- .Call("representative", .in.edges, PACKAGE = "pcalg")
          result <- new(score$.pardag.class, nodes = .nodes, in.edges = in.edges)
          result$.params <- score$global.mle(result)
          
          return(result)
        }
        ))

#' Coercion to a graphNEL instance
setAs("ess.graph", "graphNEL",
    def = function(from) {
      result <- new("graphNEL", 
          nodes = from$.nodes, 
          edgeL = from$.in.edges, 
          edgemode = "directed")
      return(reverseEdgeDirections(result))
    })

#' Coercion to a (logical) matrix
setAs("ess.graph", "matrix",
    def = function(from) {
      p <- from$node.count()
      sapply(1:p, function(i) 1:p %in% from$.in.edges[[i]])
    })

#' Plot method (needs Rgraphviz to work!!)
setMethod("plot", "ess.graph", 
    function(x, y, ...) {
      if (!validObject(x))
        stop("The parametric DAG model to be plotted is not valid")
      
      if (missing(y))
        y <- "dot"
      invisible(plot(as(x, "graphNEL"), y, ...))
    })
        
#' Gaussian causal model
setRefClass("gauss.pardag",
    contains = "pardag",
    
    validity = function(object) {
      if (any(names(object$.params) != object$.nodes))
        return("The elements of 'params' must be named after the nodes.")
      if (!all(sapply(1:object$node.count(), 
          function(i) length(object$.params[[i]]) == length(object$.in.edges[[i]]) + 2)))
        return("The number of parameters does not match the number of in-edges.")
      
      return(TRUE)
    },
    
    methods = list(
        
        #' Yields the intercept
        intercept = function() {
          sapply(.params, function(par.vec) par.vec[2])
        },
        
        #' Sets the intercept
        set.intercept = function(value) {
          for (i in 1:node.count())
            .params[[i]][2] <<- value[i]
        },
        
        #' Yields the error variances
        err.var = function() {
          sapply(.params, function(par.vec) par.vec[1])
        },
        
        #' Sets the error variances
        set.err.var = function(value) {
          for (i in 1:node.count())
            .params[[i]][1] <<- value[i]
        },
        
        #' Yields the weight matrix
        #' 
        #' TODO add a method for sparse matrices...
        weight.mat = function() {
          ## Fill in weights
          p <- node.count()
          result <- matrix(0, p, p)
          for (i in 1:p)
            result[.in.edges[[i]], i] <- .params[[i]][-c(1, 2)]
          
          ## Set row and column names
          rownames(result) <- .nodes
          colnames(result) <- .nodes
          
          return(result)
        },
        
        #' Yields an observational or interventional covariance matrix
        #' 
        #' @param   target    intervention target
        #' @param   ivent.var variances of the intervention variables
        #' @return  (observational or interventional) covariance matrix
        cov.mat = function(target = integer(0), ivent.var = numeric(0)) {
          A <- -weight.mat()
          A[, target] <- 0
          diag(A) <- 1
          A <- solve(A)
          
          all.var <- err.var()
          all.var[target] <- ivent.var
          
          return(t(A) %*% diag(all.var) %*% A)
        }
        )
    )


    