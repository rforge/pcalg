library(pcalg)
source("/u/kalischm/research/packages/pcalg/branches/dev/R/pcalg.R")
require(RBGL)

set.seed(123)

g <- randomDAG(4,0.6)
dat <- rmvDAG(1000,g,errDist="normal")

res1 <- round(pcorOrder(3,4,c(2,1),cor(dat)),8)
res2 <- round(pcorOrder(3,4,2,cor(dat)),8)

## TODO Markus: this always fails; please fix
if((res1!=-0.00474175) | (res2!=0.00890934)) {
   stop("Test of pcorOrder: Theoretical result not matched!")
}
