.onLoad <- function(libname, pkgname) {
    options(setarwav.mode="rcpp")
}

.onUnload <- function(libpath) {
    options(setarwav.mode=NULL)
}
