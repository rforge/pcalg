library(pcalg)
##################################################
## pcAlgo object
##################################################
## Load predefined data
data(gmG)
n <- nrow    (gmG8$x)
V <- colnames(gmG8$x)

## define sufficient statistics
suffStat <- list(C = cor(gmG8$x), n = n)
## estimate CPDAG
skel.fit <- skeleton(suffStat, indepTest = gaussCItest,
             alpha = 0.01, labels = V)
resSkel <- displayAmat(skel.fit)
checksSkel <- rep(FALSE,4)
checksSkel[1] <- (resSkel$type == "amat.cpdag")
checksSkel[2] <- (resSkel$amat["Author", "Bar"] == 1)
checksSkel[3] <- (resSkel$amat["Bar","Author"] == 1)
checksSkel[4] <- (resSkel$amat["Ctrl","Author"] == 0)
stopifnot(all(checksSkel))

pc.fit <- pc(suffStat, indepTest = gaussCItest,
             alpha = 0.01, labels = V)
resPC <- displayAmat(pc.fit)
checksPC <- rep(FALSE,4)
checksPC[1] <- (resPC$type == "amat.cpdag")
checksPC[2] <- (resPC$amat["V5", "V8"] == 0)
checksPC[3] <- (resPC$amat["V8","V5"] == 1)
checksPC[4] <- (resPC$amat["Goal","Author"] == 0)
stopifnot(all(checksPC))
##################################################
## fciAlgo object
##################################################
set.seed(42)
p <- 7
## generate and draw random DAG :
myDAG <- randomDAG(p, prob = 0.4)

## find PAG using the FCI algorithm
myC <- cov2cor(trueCov(myDAG))
suffStat <- list(C = myC, n = 10^9)
V <- LETTERS[1:p] ## labels of nodes

tmpFCI <- fci(suffStat, indepTest=gaussCItest, labels = V, 
           alpha = 0.9999, doPdsep = FALSE)
resFCI <- displayAmat(tmpFCI)

checksFCI <- rep(FALSE,4)
checksFCI[1] <- (resFCI$type == "amat.pag")
checksFCI[2] <- (resFCI$amat["B","E"] == 2)
checksFCI[3] <- (resFCI$amat["C","D"] == 1)
checksFCI[4] <- (resFCI$amat["G","A"] == 3)
stopifnot(all(checksFCI))
