parzen.window <- function(M, len=M) {
    if (len < M) stop("len < M")
    l <- numeric(len+1)
    r.half <- floor(M/2)
    r1 <- 0:r.half
    r2 <- (r.half+1):M
    l[r1+1] <- 1 - 6*(abs(r)/M)^2 + 6*(abs(r)/M)^3
    l[r2+1] <- 2*(1-abs(r)/M)^3
    return(l)
}
