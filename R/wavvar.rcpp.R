up.sample <- function(x, by) {
    .Call("upsample", as.numeric(x), as.integer(by))
}

convolve.rcpp <- function(x, y) {
    .Call("convolve", as.numeric(x), as.numeric(y))
}

betahat <- function(idx, Lj, j) {
    .Call("betahat", as.integer(idx), as.integer(Lj), as.integer(j))
}

modwt.missing.rcpp <- function(x, wavelet, n.levels, idx) {
    N <- length(x)
    stopifnot(N == length(idx))

    wav <- wavDaubechies(wavelet)
    h <- wav$wavelet
    g <- wav$scaling

    max.level = wavMaxLevel(length(h), N)
    stopifnot(n.levels <= max.level)

    .Call("modwt2missing", 
       as.numeric(x), h, g, as.integer(n.levels), as.integer(idx));


}
