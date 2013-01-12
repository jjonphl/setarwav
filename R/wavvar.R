# R implementation of wavelet variance computation with missing data

# utility functions
Lj <- function(L, j) (2^j-1)*(L-1)+1

# use wmtsa acvs
acvs.if <- function(x, biased=TRUE, recenter=FALSE) {
    .Call("RS_math_acvs", as.numeric(x),
                          as.logical(biased),
                          as.logical(recenter),
          PACKAGE="ifultools")
}

# compute beta
# where beta[l1,l2]^{-1} = \sum_{t=Lj-1}^{N-1} idx[t-l1]*idx[t-l2]
#  (base 0)
beta <- function(idx, Lj, j) {
    N <- length(idx)
    ret <- matrix(0, nr=Lj, nc=Lj)
    i <- 1 + seq.int(Lj,N)

    for (l1 in seq.int(Lj)) {
        for (l2 in seq.int(l1,Lj)) {
            ret[l1,l2] <- ret[l2,l1] <- sum(idx[i-l1]*idx[i-l2])
        }
    }

    zeros <- (ret == 0)
    if (any(zeros)) warning(paste("beta=0 at Lj=", Lj, ", J=", j, " consider lowering number of levels.", sep=""))
    #print(zeros)

    beta <- (N - Lj + 1) * 1/ret
    beta[zeros] <- 1             # kill those going to infinity
    beta
    #ret
}

beta.cmp <- cmpfun(beta)

# outer product of 2 vectors
# used to compte h[j,l1]*h[j,l2]
# l1,l2 = 0, ..., Lj-1
outer.product <- function(h) {
    matrix(h,nc=1) %*% matrix(h,nr=1)
}

# compute delta's:  idx[t-l1]*idx[t-l2]
# l1,l2 = 0, ..., Lj-1 (0-based)
# note that these are "reversed" sequences wrt l1,l2
delta <- function(idx, Lj, t) {
    stopifnot(t >= Lj)
    i <- idx[t - (1:Lj) + 1]
    matrix(i,nc=1) %*% matrix(i,nr=1)
}

# compute combination of X's (original data)
# to estimators are available: 
# - (X[t-l1] - X[t-l2])^2 (default)
# - X[t-l1] * X[t-l2]
# note that this sequence is also reversed, like for delta
x.op <- function(x, Lj, t, f=function(x1,x2) (x1-x2)^2) {
    stopifnot(t >= Lj)
    ret <- matrix(0, nr=Lj, nc=Lj)
    x <- x[t - (1:Lj) + 1]
    for (l1 in 1:Lj) {
        for (l2 in seq.int(l1,Lj)) {
            ret[l1,l2] <- ret[l2,l1] <- (x[l1]-x[l2])^2 #f(x[l1],x[l2])
        }
    }
    ret
}

x.op2 <- function(x, Lj, time) {
    stopifnot(time >= Lj)
    x <- x[time - (1:Lj) + 1]
    xm <- matrix(rep(x, times=Lj), nr=Lj)
    (t(xm) - xm)^2
}
x.op2.cmp <- cmpfun(x.op2)

# up-sample a vector, insert u zeros in b/w elements
up.sample <- function(x, u) {
    N <- length(x)
    ret <- rep(0, (N-1)*(u+1) + 1)
    ret[seq.int(1,length(ret), by=u+1)] <- x
    ret
}

modwt.missing.r <- function(x, wavelet, n.levels, idx) {
    N <- length(x)
    levels <- paste("d", 1:n.levels, sep="")
    filter <- wavDaubechies(wavelet, normalize=TRUE)
    h <- filter$wavelet
    g <- filter$scaling

    stopifnot(Lj(filter$length, n.levels) < N)
    ret <- list()


    for (j in seq.int(n.levels)) {
        if (j == 1) {
            hj <- h
            gj <- g
        } else {
            .g <- gj
            hj <- conv(.g, up.sample(h, 2^(j-1)-1))
            gj <- conv(.g, up.sample(g, 2^(j-1)-1))
        }
        Lj <- length(hj)

        Hj <- outer.product(hj)
        beta.j <- beta.cmp(idx, Lj)

        w <- rep(0, N-Lj+1)
        for (t in Lj:N) {
            x.j <- x.op2.cmp(x, Lj, t)
            d.j <- delta(idx, Lj, t)

            stopifnot(all(dim(x.j) == c(Lj, Lj)))
            stopifnot(all(dim(d.j) == c(Lj, Lj)))

            w[t-Lj+1] <- -0.5 * sum(Hj*beta.j*x.j*d.j)
        }

        ret[[levels[j]]] <- w

    }

    ret
}
