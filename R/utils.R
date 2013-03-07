up.sample <- function(x, by) {
    .Call("upsample", as.numeric(x), as.integer(by))
}

convolve.rcpp <- function(x, y) {
    .Call("convolve", as.numeric(x), as.numeric(y))
}


