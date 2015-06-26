#### Code for jointIda (Preetam)
#### extract.parent.sets (version Nov 14)

extract.parent.sets <- function(x.pos, amat.cpdag, isCPDAG=FALSE) {
  amat.cpdag[which(amat.cpdag != 0)] <- 1
  amat.undir <- amat.cpdag*t(amat.cpdag)
  amat.dir <- amat.cpdag - amat.undir
  pasets.dir <- lapply(x.pos,function(x) which(amat.dir[,x]!=0))

  ## get all important connected components of the undirected subgraph
  conn.comp.imp <- NULL
  x.temp <- x.pos
  while (length(x.temp)>0){
    comp.temp <- dfs(graph = graph.adjacency(amat.undir,mode='undirected'),
                     root = x.temp[1], unreachable=FALSE)$order
    comp.temp <- comp.temp[!is.na(comp.temp)]
    x.temp <- setdiff(x.temp,comp.temp)
    conn.comp.imp <- c(conn.comp.imp,list(comp.temp))
  }
  ## new igraph 1.0.0 has more than just node numbers.  need only those
  ## "Hack":
  conn.comp.imp <- lapply(conn.comp.imp, as.integer)

  ## Chordality test, if required
  chordal <- if (!isCPDAG) {
    chrordality.test.fun <- function(comp)
      is.chordal(graph.adjacency(amat.undir[comp,comp],mode="undirected"),fillin=F)$chordal
    sapply(conn.comp.imp, chrordality.test.fun)
  }else{
    rep(TRUE,length(conn.comp.imp))
  }

  ## Function for getting locally valid parent sets
  all.locally.valid.parents.undir <- function(amat,x){ # x must be a scaler
    pa.dir <- pasets.dir[[which(x.pos==as.integer(rownames(amat)[x]))]]
    paset <- list(pa.dir)
    pa <- which(amat[,x] != 0) # cannot be a null set
    if (length(pa)==1){
      paset <- c(paset,list(c(as.integer(rownames(amat)[pa]),pa.dir)))
    }else{
      for (i in 1:length(pa)){
        pa.tmp <- combn(pa, i, simplify = TRUE)
        n.comb <- ncol(pa.tmp)
        for (j in 1:n.comb) {
          pa.t <- pa.tmp[, j]
          new.coll <- FALSE
          if (length(pa.t)>1){
            tmp <- amat[pa.t,pa.t]
            diag(tmp) <- 1
            new.coll <- (min(tmp)==0)
          }
          if (!new.coll) paset <- c(paset,list(c(as.integer(rownames(amat)[pa.t]),pa.dir)))
        }
      }
    }
    return(paset)
  }


  extract.parent.sets.from.conn.comp <- function(i){
    all.nodes <- conn.comp.imp[[i]]
    nvar <- length(all.nodes)
    if (nvar==1){
      pasets.comp <- list(pasets.dir[match(all.nodes,x.pos)])
    }else{
      conn.comp.mat <- amat.undir[all.nodes,all.nodes]
      rownames(conn.comp.mat) <- all.nodes
      x.temp <- intersect(all.nodes,x.pos)
      x.temp2 <- match(x.temp,all.nodes)
      if(chordal[i] & nvar<=12){
        rownames(conn.comp.mat) <- colnames(conn.comp.mat) <- 1:nvar
        all.extensions <- allDags(conn.comp.mat,conn.comp.mat,NULL)
        pa.fun <- function(amat,j) c(all.nodes[which(amat[,x.temp2[j]]!=0)],
                                     pasets.dir[[match(x.temp[j],x.pos)]])
        parent.sets.fun <- function(r) lapply(1:length(x.temp), pa.fun,
                                              amat=matrix(all.extensions[r,],nrow=nvar))
        pasets.comp <- lapply(1:nrow(all.extensions),parent.sets.fun)
      }else{
        pasets.comp <- lapply(x.temp2,all.locally.valid.parents.undir,amat=conn.comp.mat)
        idx <- expand.grid(lapply(1:length(x.temp),function(j) 1:length(pasets.comp[[j]])))
        pasets.comp <- lapply(1:nrow(idx),function(r)
          lapply(1:length(x.temp), function(j) pasets.comp[[j]][[idx[r,j]]]))
      }

    }
    return(pasets.comp)
  }

  all.pasets <- lapply(1:length(conn.comp.imp),extract.parent.sets.from.conn.comp)
  idx <- expand.grid(lapply(1:length(all.pasets),function(i) 1:length(all.pasets[[i]])))
  x.conn.comp <- unlist(lapply(1:length(conn.comp.imp),function(i) intersect(conn.comp.imp[[i]],x.pos)))
  all.pasets <- lapply(1:nrow(idx),function(i) unlist(lapply(1:length(conn.comp.imp),
                                                             function(j) all.pasets[[j]][[idx[i,j]]]),recursive=F)[match(x.pos,x.conn.comp)])

  return(all.pasets)
}
## jointIda
jointIda <- function(x.pos,y.pos,mcov,graphEst=NULL,all.pasets=NULL,technique=c("RRC","MCD")){

  if (is.null(all.pasets)){
    amat <- as(graphEst,"matrix")
    amat[which(amat != 0)] <- 1
    all.pasets <- extract.parent.sets(x.pos,amat)
  }else{
    correct.format <- TRUE
    if (!is.list(all.pasets)) correct.format <- FALSE
    for (i in 1:length(all.pasets)){
      if (length(all.pasets[[i]]) != length(x.pos)) correct.format <- FALSE
    }
      if (!correct.format) stop("all.pasets is not given in an appropriate format.")
  }
  if (length(y.pos)>1){
    return(lapply(y.pos,function(y) jointIda(x.pos,y,mcov,all.pasets=all.pasets,technique=technique)))
  }else{
    technique <- match.arg(technique)
    if (is.element(y.pos,x.pos)) return(matrix(0,nrow=length(x.pos),ncol=length(all.pasets)))
    joint.effects <- switch(technique,
                            RRC = matrix(unlist(lapply(all.pasets,function(pasets) RRC(mcov,x.pos,y.pos,pasets))),
                                        nrow=length(x.pos)),
                            MCD = matrix(unlist(lapply(all.pasets,function(pasets) MCD(mcov,x.pos,y.pos,pasets))),
                                         nrow=length(x.pos)))
    return(joint.effects)
  }
}

## MCD
MCD <- function(cov.mat,intervention.set,var.y,pasets,return.modified.cov.mat=FALSE) {
  if (is.element(var.y,intervention.set) & !return.modified.cov.mat)
    return(rep(0,length(intervention.set)))
  if (length(intervention.set)==1 & is.element(var.y,unlist(pasets)) & !return.modified.cov.mat) return(0)

  if (!return.modified.cov.mat){
    imp.var <- unique(c(intervention.set,unlist(pasets),var.y))
    if (length(imp.var) < nrow(cov.mat)){
      cov.mat <- cov.mat[imp.var,imp.var]
      intervention.set <- match(intervention.set,imp.var)
      var.y <- match(var.y,imp.var)
      pasets <- if(length(intervention.set)>1) lapply(pasets,function(x) match(unlist(x),imp.var)) else match(unlist(pasets),imp.var)
    }
  }

  do.Cholesky.modification <- function(x){
    pa.x <-  if(length(intervention.set)>1) pasets[[match(x,intervention.set)]] else unlist(pasets)
    if (length(pa.x)==0) return(cov.mat)
    ind <- c(pa.x,x,(1:nrow(cov.mat))[-c(pa.x,x)])
    cov.mat <- cov.mat[ind,ind]
    x <- match(x,ind)
    pa.x <- match(pa.x,ind)
    temp <- gchol(cov.mat)
    Lower.tri.mat <- solve(as.matrix(temp))
    tmp1 <- bdsmatrix::diag(temp)
    Diag.mat <- base::diag(tmp1)
    Lower.tri.mat[x,pa.x] <- 0
    cov.mat <- solve(Lower.tri.mat)%*%Diag.mat%*%t(solve(Lower.tri.mat))
    return(cov.mat[order(ind),order(ind)])
  }
  for (i in 1:length(intervention.set)) {
    cov.mat <- do.Cholesky.modification(intervention.set[i])
  }

  if(return.modified.cov.mat) return(cov.mat)
  MCD.estimate <- function(x) {
    if (is.element(var.y,unlist(pasets[match(x,intervention.set)]))) {
      return(0)
    } else {
      return(cov.mat[var.y,x]/cov.mat[x,x])
    }
  }

  return(sapply(intervention.set,MCD.estimate))

}

## RRC
RRC <- function(cov.mat,intervention.set,var.y,pasets) {

  adjusted.regression <- function(x,y){
    if (x == y) return(0)
    pa.x <-  if(length(intervention.set)>1) pasets[[match(x,intervention.set)]] else unlist(pasets)
    if(is.element(y,pa.x)) return(0)
    return(solve(cov.mat[c(x,pa.x),c(x,pa.x)],cov.mat[c(x,pa.x),y])[1])
  }

  ## Define the vector of causal effects of intervention variables on var.y
  intervention.set.on.var.y <- sapply(intervention.set,adjusted.regression,y=var.y)

  ## Define the required matrix of single intervention effects
  if (length(intervention.set)>1){
    intervention.set.on.intervention.set <- matrix(apply(expand.grid(intervention.set,intervention.set),
                                                         1,function(x) adjusted.regression(x[1],x[2])),nrow=length(intervention.set))
  }else{
    return(intervention.set.on.var.y)
  }

  joint.effect.fun  <- function(x){
    if(is.element(var.y,unlist(pasets[match(x,intervention.set)]))) return(0)
    x.temp  <- match(x,intervention.set)
    intervention.set.temp  <- intervention.set[-x.temp]
    ## Intitialize the RR estimate as the single intervention effect of intervention.set on var.y
    RR.estimate <-  intervention.set.on.var.y[x.temp]
    ## Define the vector of causal effects of intervention.set on other intervention variables
    x.temp.on.intervention.set.temp <- intervention.set.on.intervention.set[x.temp,-x.temp]
    ## Define the vector of causal effects of "other" intervention variables on var.y
    intervention.set.on.var.y.temp <- intervention.set.on.var.y[-x.temp]
    ## Define the required matrix of single intervention effects
    intervention.set.on.intervention.set.temp <- intervention.set.on.intervention.set[-x.temp,-x.temp]

    while (length(intervention.set.temp)>1){
      ## update RR.estimate and the other things accounting for the elimination of the first entry of the current intervention set
      RR.estimate <- RR.estimate - x.temp.on.intervention.set.temp[1]*intervention.set.on.var.y.temp[1]
      intervention.set.on.var.y.temp <- intervention.set.on.var.y.temp[-1] - intervention.set.on.intervention.set.temp[-1,1] * intervention.set.on.var.y.temp[1]
      x.temp.on.intervention.set.temp <- x.temp.on.intervention.set.temp[-1] - x.temp.on.intervention.set.temp[1]*intervention.set.on.intervention.set.temp[1,-1]
      if (length(intervention.set.temp)>2)
        intervention.set.on.intervention.set.temp <- intervention.set.on.intervention.set.temp[-1,-1] - matrix(intervention.set.on.intervention.set.temp[-1,1],ncol=1)%*%matrix(intervention.set.on.intervention.set.temp[1,-1],nrow=1)
      intervention.set.temp <- intervention.set.temp[-1]
    }
    return(RR.estimate-x.temp.on.intervention.set.temp*intervention.set.on.var.y.temp)
  }

  return(sapply(intervention.set,joint.effect.fun))
}


###  MM: (ess-set-style 'DEFAULT) : we have much nesting ==> only indent by 2
## Local Variables:
## eval: (ess-set-style 'DEFAULT 'quiet)
## delete-old-versions: never
## End:

