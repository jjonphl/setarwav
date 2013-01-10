bp.window <- function(M, len=M) {
    if (len < M) stop("len < M")
    l <- numeric(len+1)
    r <- 1:M   # no 0
    l[1] <- 1
    l[r+1] <- 3*M^2/(pi*r)^2 * (sin(pi*r/M) / (pi*r/M) - cos(pi*r/M))
    return(l)
}
