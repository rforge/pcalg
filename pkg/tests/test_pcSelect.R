library(pcalg)

(corMeths <- eval(formals(mcor)$method))
## the other methods (but the default):
(cMeths <- corMeths[-1])

p <- 10
n <- 1000

set.seed(101)
myDAG <- randomDAG(p, prob = 0.2)
d.mat <- rmvDAG(n, myDAG, errDist = "normal")
y <- d.mat[,10]
dm <- d.mat[,-10]
res1 <- pcSelect(y,dm, alpha=0.05)
if (!all(res1$G == 1:9 %in% c(4,5,6)))
  stop("Test of pcSelect: Consistency problem 101")
dput(signif(res1$zMin, 7))
zM.1 <- c(0.07643311, 0.6568222, 0.7363025,
          11.62001, 6.621493, 18.64382,
          1.447, 1.212824, 1.221812)
stopifnot(all.equal(res1$zMin, zM.1, tol = 1e-6))

## Now all other methods:
zMin. <- setNames(as.list(cMeths), cMeths)
for(meth in cMeths) {
    cat(meth,":\n")
    rr <- pcSelect(y,dm, alpha=0.05, corMethod = meth)
    iG <- which(as.vector(rr$G))
    if (!identical(iG, c(4L,5L,6L))) 
        stop(sprintf("pcSelect(dm101.., corMethod=\"%s\"): which(G) = %s\n",
                     meth, paste(iG, collapse=", ")))
    zMin.[[meth]] <- rr$zMin
}

## dput(lapply(zMin., signif, digits=7))
zMin1.Exp <- list(
    Qn = c(0.5316588, 1.088152, 0.8118473, 11.34048, 
           6.416001, 18.79211, 0.9179979, 0.4779155, 1.346989),
    QnStable = c(0.9129029, 1.710694, 0.8471604, 10.56123, 6.168661,
                 18.88213, 0.9710747, 0.4337806, 1.140709),
    ogkScaleTau2 = c(0.4698467, 0.886344, 0.5731621, 10.94707, 6.335696,
                     18.60215, 0.9969964, 1.045772, 1.244915), 
    ogkQn = c(0.7588748, 0.6270149, 0.6489765, 11.19786, 6.052107, 
              18.683, 0.8969601, 0.9808662, 1.118423),
    shrink = c(0.07570146, 0.6493067, 0.7303193, 11.49862, 6.556261,
               18.45523, 1.433129, 1.199766, 1.211877))
stopifnot(all.equal(zMin., zMin1.Exp, tol = 1e-6))

set.seed(456)
myDAG <- randomDAG(p, prob = 0.2)
d.mat <- rmvDAG(n, myDAG, errDist = "normal")
y <- d.mat[,10]
dm <- d.mat[,-10]
iG.Exp <- c(5L,8L,9L)
res2 <- pcSelect(y,dm,alpha=0.05)
if (!identical(which(as.vector(res2$G)), iG.Exp))
  stop("Test of pcSelect: Consistency problem 456")
## dput(signif(res2$zMin, 7))
zM.2 <- c(0.1469954, 0.915466, 0.5442373,
          0.1089084, 14.16709, 1.056202, 
          0.4840786, 17.10435, 14.51551)
stopifnot(all.equal(res2$zMin, zM.2, tol = 1e-6))

## Now all other methods:
for(meth in cMeths) {
    cat(meth,":\n")
    rr <- pcSelect(y,dm, alpha=0.05, corMethod = meth)
    iG <- which(as.vector(rr$G))
    if (!identical(iG, iG.Exp))
        (if(all(iG.Exp %in% iG)) message else stop)(
            sprintf("pcSelect(dm456.., corMethod=\"%s\"): which(G) = %s\n",
                    meth, paste(iG, collapse=", ")))
    zMin.[[meth]] <- rr$zMin
}

## dput(lapply(zMin., signif, digits=7))
zMin2.Exp <- 
    list(
        Qn = c(1.055523, 2.555242, 1.691455, 0.7796301, 11.39895,
               0.8692785, 0.04107331, 14.0156, 14.82198),
        QnStable = c(0.6208816, 2.475236, 1.831045, 0.9328639, 11.11562,
                     1.139745, 0.04850462, 14.21195, 14.77496),
        ogkScaleTau2 = c(0.4164263, 1.092748, 1.689739, 0.4691277, 11.92243,
                         1.028907, 0.06462646, 13.80103, 14.60776),
        ogkQn = c(0.237811, 1.916074, 1.704892, 0.3637768, 14.42101, 
                  0.8777305, 0.05152724, 16.72384, 14.60674),
        shrink = c(0.1459792, 0.902588, 0.540349, 0.1081555, 13.98006,
                   1.047259, 0.4807316, 16.88456, 14.4006))
stopifnot(all.equal(zMin., zMin2.Exp, tol = 1e-6))
