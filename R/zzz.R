.First.lib <- function(libname, pkgname) {
    options(setarwav=list(modwt.missing="r"))
}

.Last.lib <- function() {
    options(setarwav=NULL)
}
