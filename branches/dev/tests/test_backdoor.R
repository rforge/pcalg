library(pcalg)


#####################################################################
##CPDAG
#####################################################################
##################################################
## Example not identifiable
## Maathuis and Colombo (2013), Fig. 3, p.14
##################################################

## create the graph
p <- 5
amat <- t(matrix(c(0,0,1,1,1, 0,0,1,1,1, 0,0,0,1,0, 0,0,0,0,1, 0,0,0,0,0),5,5))
colnames(amat) <- rownames(amat) <- as.character(1:5)
V <- as.character(1:5)
edL <- vector("list",length=5)
names(edL) <- V
edL[[1]] <- list(edges=c(3,4,5),weights=c(1,1,1))
edL[[2]] <- list(edges=c(3,4,5),weights=c(1,1,1))
edL[[3]] <- list(edges=4,weights=c(1))
edL[[4]] <- list(edges=5,weights=c(1))
g <- new("graphNEL", nodes=V, edgeL=edL, edgemode="directed")

## estimate the true CPDAG
myCPDAG <- dag2cpdag(g)
## Extract the adjacency matrix of the true CPDAG
true.amat <- as(myCPDAG, "matrix")
true.amat[which(true.amat!=0)] <- 1

## The effect is not identifiable, in fact:
tmp.set <- backdoor(true.amat, 3, 5, type="cpdag")

if (!is.na(tmp.set)) {
  stop("Test of backdoor: CPDAG example wrong set found!")
}


##################################################################
##PAG
##################################################################
##################################################
## Example identifiable
## Maathuis and Colombo (2013), Fig. 7, p.17
##################################################

## create the graph
p <- 7
amat <- t(matrix(c(0,0,1,1,0,0,0, 0,0,1,1,0,0,0, 0,0,0,1,0,1,0, 0,0,0,0,0,0,1, 0,0,0,0,0,1,1, 0,0,0,0,0,0,0, 0,0,0,0,0,0,0),7,7))
colnames(amat) <- rownames(amat) <- as.character(1:7)
V <- as.character(1:7)
edL <- vector("list",length=7)
names(edL) <- V
edL[[1]] <- list(edges=c(3,4),weights=c(1,1))
edL[[2]] <- list(edges=c(3,4),weights=c(1,1))
edL[[3]] <- list(edges=c(4,6),weights=c(1,1))
edL[[4]] <- list(edges=7,weights=c(1))
edL[[5]] <- list(edges=c(6,7),weights=c(1,1))
g <- new("graphNEL", nodes=V, edgeL=edL, edgemode="directed")
L <- 5

## compute the true covariance matrix of g
cov.mat <- trueCov(g)

## transform covariance matrix into a correlation matrix
true.corr <- cov2cor(cov.mat)
suffStat <- list(C=true.corr, n=10^9)
indepTest <- gaussCItest

## estimate the true PAG
true.pag <- dag2pag(suffStat, indepTest, g, L, alpha = 0.9999)

## The effect is identifiable and
tmp.set <- backdoor(true.pag@amat, 4, 6, type="pag")
true.set <- c(1,2)

if (!all(tmp.set==true.set)) {
  stop("Test of backdoor: PAG example wrong set found!")
}
