library(pcalg)

m <- matrix(FALSE, 3,3)
tmp <- list(NULL, NULL, NULL)
sepset <- list(tmp, tmp, tmp)
mm <- mode(udag2pag(m,sepset))
if (mm != "numeric") {
    stop("Test of udag2pag: Output does not have mode numeric!")
}
