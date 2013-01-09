setar.sim.simple <-
function(thresholds, model, n, delay, 
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
