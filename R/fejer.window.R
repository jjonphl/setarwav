fejer.window <- function(M, len=M) {
    if (len < M) stop("len < M")
    l <- numeric(len+1)
    r <- 0:M
    l[r+1] <- 1 - abs(r)/M
    return(l)
}
