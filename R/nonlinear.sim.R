# non-linear simulators

# SETAR - does not compute "debugging info" like regime per obs
# general form: 
# model=list()
setar.sim.simple <- function(n, thresholds, model, delay, 
                             sd=1, startup=100){
    if (!is.list(model)) {
        stop("model should be a list(ar.regime1, ar.regime2, ...)")
    }

    if (length(thresholds)+1 != length(model)) {
        stop("Length of thresholds & model params does not match.")
    }

    if (any(diff(thresholds) < 0)) {
        stop("Thresholds must be specified in ascending order")
    }

    len <- startup + n      # length including burn-in
    rs <- rep(0, len)       # to record regime per time, 1-based
    thresholds <- c(-Inf, sort(thresholds), Inf)
    n.regimes <- length(thresholds) + 1
    if (length(sd) != n.regimes) sd <- rep(sd, length=n.regimes)

    y <- numeric(len)       # pre-allocate setar sim output
    e.list <- list()
    for (i in 1:n.regimes) {
        # for cases where different regimes have different error sd
        e.list[[i]] <- rnorm(len, sd=sd[i])
    }

    for (i in 2:len) {
        d <- max(i-delay, 1)

        # compute setar
        regime <- which(y[d] < thresholds)[1] - 1 
        coefs <- model[[regime]]
        p <- length(coefs) - 1         # 1st arg is intercept
        idx <- seq(i-1, max(i-p,1))
        lags <- c(1, y[idx], rep(0,p-length(idx)))
        y[i] <- sum(lags*coefs) + e.list[[regime]][i]

    }

    return(y[(startup+1):len])
}

# same name as tsDyn::setar.sim
# PARAMS
#   thresholds - array of threshold values
#   model - list of AR coefficients (p0, p1, ..., pp), p0 = intercept
#   n - number of obs to produce
#   delay - delay param
#   sd - stdev of the generated innovations. sd=0 gives the skeleton
#   n.start - 
#
setar.sim2 <- function(thresholds, model, n, delay, 
                      sd=1.5*(max(thresholds)-min(thresholds)), 
                      startup=100, as.dataframe=FALSE){
    if (!is.list(model)) {
        stop("model should be a list(ar.regime1, ar.regime2, ...)")
    }

    if (length(thresholds)+1 != length(model)) {
        stop("Length of thresholds & model params does not match.")
    }

    if (any(diff(thresholds) < 0)) {
        stop("Thresholds must be specified in ascending order")
    }

    len <- startup + n      # length including burn-in
    skel <- rep(0, len)     # "skeleton" (not strictly), without error term
    rs <- rep(0, len)       # to record regime per time, 1-based
    thresholds <- c(-Inf, sort(thresholds), Inf)
    n.regimes <- length(thresholds) + 1
    if (length(sd) != n.regimes) sd <- rep(sd, length=n.regimes)

    y <- numeric(len)       # will hold simulated setar realization
    e.list <- list()
    for (i in 1:n.regimes) {
        # for cases where different regimes have different error sd
        e.list[[i]] <- rnorm(len, sd=sd[i])
    }

    #if (all(sd == 0)) {  
    #    # if sd == 0, caller is requesting skeleton only
    #    default.sd <- 1.5*(max(thresholds)-min(thresholds))
    #    y[1:startup] <- rnorm(startup, mean(thresholds), default.sd)
    #}

    for (i in 2:len) {
        d <- max(i-delay, 1)

        # compute skeleton - skel regime may be diff from setar regime below!
        regime <- which(skel[d] < thresholds)[1] - 1
        coefs <- model[[regime]]
        p <- length(coefs) - 1    # coefs include intercept
        idx <- seq(i-1, max(i-p,1))
        lags <- c(1, skel[idx], rep(0, p-length(idx)))
        skel[i] <- sum(lags*coefs) + (i <= startup) * y[i]

        # compute setar
        regime <- which(y[d] < thresholds)[1] - 1 
        rs[i] <- regime              # record previous regime
        coefs <- model[[regime]]
        p <- length(coefs) - 1         # 1st arg is intercept
        idx <- seq(i-1, max(i-p,1))
        lags <- c(1, y[idx], rep(0,p-length(idx)))
        y[i] <- sum(lags*coefs) + e.list[[regime]][i]

    }

    # record last regime after exit
    if (as.dataframe) {
        regime[n] <- which(y[d] < thresholds)[1] - 1 

        ret.df <- data.frame(regime=as.integer(rs[(startup+1):len]), 
                             skel=skel[(startup+1):len],
                             y=y[(startup+1):len]);
        attr(ret.df, "call") <- sys.call()
        return(ret.df)
    } else {
        list(y=ts(y[(startup+1):len]), skel=ts(skel[(startup+1):len]),
             regimes=as.integer(rs[(startup+1):len]), call=sys.call())
    }
}

# EXPAR
# general form:
# X[t] = \sum_{i=1}^p (a[i] + (b[i] + c[i]*X[t-d])*exp(-d[i]*X[t-d]^2)) *
#                      X[t-i] + e[t]
# model = list(a=(a1,a2,...,a[p]),b=(b1,b2,...,b[p]),
#              c=(c1,c2,...,c[p]),d=(d1,d2,...,d[p]))
expar.sim <- function(n, model, d, sd=1, startup=100) {
    if (!is.list(model) || !setequal(names(model), c("a","b","c","d"))) {
        stop("model should be a list(a=...,b=...,c=...,d=...)")
    }

    pr <- sapply(model, length)
    p <- min(pr)
    if (p != max(pr)) {
        stop("all model vectors should have the same length")
    }

    len <- startup + n      # length including burn-in
    y   <- numeric(n)
    e   <- rnorm(len, sd=sd)

    for (i in 2:len) {
        d0 <- max(i-d, 1)
        yd <- y[max(i-d,1)]
        ar.coefs <- with(model, a + (b + c*yd)*exp(-d*yd^2))
        idx <- seq.int(i-1, max(i-p,1))
        lags <- c(y[idx], rep(0,p-length(idx))) 
        y[i] <- sum(ar.coefs*lags) + e[i]
    }

    return(y[(startup+1):len])
}


# BILINEAR
# general form:
# X[t] = a0 + \sum_{i=1}^p a[i]*X[t-i] + e[t] + \sum_{i=1}^q b[i]*e[t-i] +
#        \sum_{i=1}^P \sum_{j=1}^Q c[i,j]*X[t-i]*e[t-j]
# model = list(a0=a0,ar=c(a1,...,a[p]), ma=(b1,...,b[q]), 
#              bl=matrix(c(k[1,1],k[2,1],...),nr=P,nc=Q))
bilinear.sim <- function(n, model, sd=1, startup=100) {
    if (!is.list(model) || any(!(names(model) %in% c("ar","ma","bl","a0")))) {
        stop("model should be a list(a0=...,ar=...,ma=...,bl=...)")
    }

    p <- length(model$ar)
    q <- length(model$ma)
    a0 <- if (!is.null(model$a0)) model$a0 else 0

    if (is.null(model$bl)) {
        warning("No BL model, might as well use arima.sim!")
        P <- Q <- 0
    } else {
        bl.mat <- as.matrix(model$bl)
        P <- nrow(bl.mat)
        Q <- ncol(bl.mat)
    }

    len <- startup + n
    y   <- numeric(n)
    e   <- rnorm(len, sd=sd)

    for (i in 2:len) {
        ar.sum <- if (p > 0) {
            idx <- seq.int(i-1, max(i-p,1))
            ar.lag <- c(y[idx], rep(0,p-length(idx)))
            sum(model$ar*ar.lag)
        } else 0
        
        ma.sum <- if (q > 0) {
            idx <- seq.int(i-1, max(i-q,1))
            ma.lag <- c(e[idx], rep(0,p-length(idx)))
            sum(model$ma*ma.lag)
        } else 0

        bl.sum <- if (P*Q > 0) {
            idx <- seq.int(i-1, max(i-P,1))
            bl.ylag <- c(y[idx], rep(0,P-length(idx)))
            idx <- seq.int(i-1, max(i-Q,1))
            bl.elag <- c(e[idx], rep(0,Q-length(idx)))
            (bl.ylag %*% bl.mat %*% bl.elag)
        } else 0

        y[i] <- a0 + ar.sum + e[i] + ma.sum + bl.sum
    }

    return(y[(startup+1):len])
}
