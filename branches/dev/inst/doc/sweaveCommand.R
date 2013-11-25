## A version that works for "everyone" (not just Markus)
## *after* having installed or already reqquire()d  'pcalg':

manualInst.vignette <- function(fstem) {
    Rnw <- system.file("doc", paste(fstem,"Rnw", sep="."), package="pcalg")
    if(interactive())
        file.show(Rnw)
    Sweave (Rnw)
    tools::texi2pdf(paste(fstem, "tex", sep="."))
    Stangle(Rnw)
    ## and test if the code wokrs:
    Rfile <- paste(fstem, "R", sep=".")
    source(Rfile)

    pkgSrcDir <-
        switch(Sys.getenv("USER"),
               "maechler" = "~/R/Pkgs/pcalg-dev",
               "kalischm" = ".....",    # PATH of pcalg-dev
               stop("Must add your username in doc/sweaveCommand.R "))
    ## Now manually "install" the vignette pdf and R :

    (pkgDoc <- file.path(pkgSrcDir, "inst","doc"))
    file.copy(c(Rnw, Rfile, paste(fstem, "pdf", sep=".")),
              pkgDoc)
    ## --> TRUE TRUE TRUE  if it works
}

stopifnot(manualInst.vignette(fstem = "pcalgDoc"))
