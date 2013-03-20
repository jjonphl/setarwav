partition <- function(y, delay, q) {
    cutoff <- quantile(y, prob=c(q, 1-q), type=4)
    nobs <- length(y)
    eff.nobs <- nobs - delay

    idx.lo <- (y[1:eff.nobs] < cutoff[2])
    idx.up <- (y[1:eff.nobs] > cutoff[1])
    #idx.up <- (y[1:eff.nobs] > cutoff[1])

    y.eff <- y[(delay+1):nobs]
    y.lo <- y.eff * as.numeric(idx.lo)
    y.up <- y.eff * as.numeric(idx.up)

    list(y.eff=y.eff,idx.lo=idx.lo,idx.up=idx.up,y.lo=y.lo,y.up=y.up)
}

wav.test <- function(y, delay, q=0.333, wavelet="s8", fill.mean=F, 
                      pvalue.alpha=0.05, stat.type=-1, 
                      modwt.missing=modwt.missing.r, #TODO: modwt.missing.rcpp
                      debug=TRUE, levels=NULL) {
    cutoff <- quantile(y, prob=c(q, 1-q), type=4)

    nobs <- length(y)
    eff.nobs <- nobs - delay

    idx.lo <- (y[1:eff.nobs] < cutoff[2])
    idx.up <- (y[1:eff.nobs] > cutoff[1])
    #idx.up <- (y[1:eff.nobs] > cutoff[1])

    y.eff <- y[(delay+1):nobs]
    y.lo <- y.eff * as.numeric(idx.lo)
    y.up <- y.eff * as.numeric(idx.up)

    stopifnot(length(y.lo) == eff.nobs)

    # center at 0
    #y.lo[idx.lo] <- y.lo[idx.lo] - sum(y.lo)/sum(idx.lo)
    #y.up[idx.up] <- y.up[idx.up] - sum(y.up)/sum(idx.up)

    # compute modwt
    wav    <- wavDaubechies(wavelet, normalize=TRUE) # MODWT is normalized
    if (is.null(levels)) {
        # do not include last level to be sure beta^{-1} does not
        # have 0 element
        levels <- wavMaxLevel(wav$length, length(y.eff), "modwt") - 1
    }

    modwt.y <- wavMODWT(y.eff, wavelet, levels)
    modwt.lo <- modwt.missing(y.lo, wavelet, levels, idx.lo)
    modwt.up <- modwt.missing(y.up, wavelet, levels, idx.up)

    # get interior of MODWT, i.e. those coefs that do not depend on boundary
    # assumptions
    int.y <- wavBoundary(modwt.y)

    # compute unbiased variance
    var.d <- function(x) sum(x^2)/length(x)

    wavcoefs <- names(modwt.y$data)[grep("^d",names(modwt.y$data))]
    var.lo <- sapply(modwt.lo, mean)  # modwt.missing output is |W[j,t]|^2
    var.up <- sapply(modwt.up, mean)

    # lengths
    Mj <- eff.nobs - Lj(wav$length, seq_along(wavcoefs)) + 1
    names(Mj) <- wavcoefs

    # "spectral energy" at each level j
    acvs.orig <- function(x, lag) {   # biased autocovariance function
        s <- rep(0, lag+1)
        N <- length(x)
        for (tau in 0:lag) {
            x1 <- x[1:(N-tau)]
            x2 <- x[(1+tau):N]
            s[tau+1] <- sum(x1*x2) / (lag+1)  # R is 1-based
        }
        return(s)
    }

    acvs <- function(x, lag) {
        # use DFT-based acvs
        acvs.if(x[1:(lag+1)], biased=TRUE, recenter=FALSE)   
    }

    Aj.est <- function(dwt, wavcoefs, Mj) {
        Aj <- rep(0, length(wavcoefs))
        for (j in seq_along(wavcoefs)) {
            x <- dwt[[wavcoefs[j]]]
            s <- acvs(x, Mj[j]-1)
            Aj[j] <- s[1]^2/2 + sum(s[2:Mj[j]]^2)
        }
        names(Aj) <- wavcoefs
        Aj
    }

    #Aj <- Aj.est(int.y$interior, wavcoefs, Mj)
    #Aj.up <- Aj.est(modwt.up, wavcoefs, Mj)
    #Aj.lo <- Aj.est(modwt.lo, wavcoefs, Mj)

    # location stat, not so good
    location.stat <- function(n) {
        res <- (var.up[n] - var.lo[n]) / sqrt((2*Aj[n]/Mj[n]))
        names(res) <- NULL
        res
    }

    # location stat, adjust with corresponding Aj est
    location.stat2 <- function(n) {
        res <- var.up[n]/sqrt(2*Aj.up[n]/Mj[n]) -
               var.lo[n]/sqrt(2*Aj.lo[n]/Mj[n])
        names(res) <- NULL
        res
    }

    # scale (ratio) stat
    scale.stat <- function(n) {
        res <- (var.up[n] / sqrt(2*Aj.up[n]/Mj[n])) /
               (var.lo[n] / sqrt(2*Aj.lo[n]/Mj[n]))
        names(res) <- NULL
        res
    }

    # to make it a 1-sided test, make sure the larger coef is at the bottom
    scale.stat2 <- function(n) {
        up <- (var.up[n] / sqrt(2*Aj.up[n]/Mj[n]))
        lo <- (var.lo[n] / sqrt(2*Aj.lo[n]/Mj[n]))

        res <- max(up,lo)/min(up,lo)
        names(res) <- NULL
        res
    }

    # Aj.up/lo is messing up the F, use common Aj (w/c cancels out)
    scale.stat3 <- function(n) {
        up <- var.up[n]
        lo <- var.lo[n]

        res <- max(up,lo)/min(up,lo)
        names(res) <- NULL
        res
    }

    sum.stat <- function() {
        up <- sum(var.up)
        lo <- sum(var.lo)

        res <- max(up,lo)/min(up,lo)
        names(res) <- NULL
        res
    }

    ratio.sum.stat <- function() {
        res <- sum(log(var.up/var.lo))
        if (is.nan(res)) print(paste(var.up, var.lo))
        names(res) <- NULL
        res
    }

    ratio.sum.stat2 <- function() {
        res <- sum(abs(log(var.up/var.lo)))
        if (is.nan(res)) print(paste(var.up, var.lo))
        names(res) <- NULL
        res
    }

    ratio.sum.stat3 <- function() {
        tmp <- abs(log(var.up/var.lo))
        if (any(is.nan(tmp))) {
            print(paste(var.up, var.lo))
            return(NA)
        }
        tmp <- tmp * 1/2^(seq.int(0,length(tmp)-1))
        res <- sum(tmp)
        names(res) <- NULL
        res
    }

    df.up <- sum(idx.up)
    df.lo <- sum(idx.lo)

    # compute p-value for location.stat2
    location.pvalue2 <- function(n) {
        res <- pnorm(abs(stat[n]), 0, 2)
        names(res) <- NULL
        res
    }

    # compte p-value for scale.stat2
    scale.pvalue2 <- function(n) {
        #up <- (var.up[n] / sqrt(2*Aj.up[n]/Mj[n]))
        #lo <- (var.lo[n] / sqrt(2*Aj.lo[n]/Mj[n]))
        up <- var.up[n]
        lo <- var.lo[n]


        res <- if (up > lo) {
            1 - pf(stat[n], df.up, df.lo)
        } else {
            1 - pf(stat[n], df.lo, df.up)
        }
        names(res) <- NULL
        res
    }

    scale.pvalue.edof1 <- function(n) {
        up <- var.up[n]
        lo <- var.lo[n]

        edof1.up <- up^2 * Mj[n] / Aj.up[n]
        edof1.lo <- lo^2 * Mj[n] / Aj.up[n]

        res <- if (up > lo) {
            1 - pf(stat[n], edof1.up, edof1.lo)
        } else {
            1 - pf(stat[n], edof1.lo, edof1.up)
        }
        names(res) <- NULL
        res
    }
    
    sum.pvalue <- function() {
        up <- sum(var.up)
        lo <- sum(var.lo)

        res <- if (up > lo) {
            1 - pf(stat, df.up, df.lo)
        } else {
            1 - pf(stat, df.lo, df.up)
        }
    }

# test notes:
# scale.stat2 + scale.pvalue2 = accepts everything (due to Aj??)
# scale.stat3 + scale.pvalue2 = rejects almost everything (no AJ)
# location.stat2 + location.pvalue2 = accepts everything

    if (stat.type == 1) {
        stat.fun <- scale.stat2
        pvalue.fun <- scale.pvalue2
        stat.list <- TRUE
    } else if (stat.type == 2) {
        stat.fun <- location.stat2
        pvalue.fun <- location.pvalue2
        stat.list <- TRUE
    } else if (stat.type == 3) {
        stat.fun <- sum.stat
        pvalue.fun <- sum.pvalue
        stat.list <- FALSE
    } else if (stat.type == 4) {
        stat.fun <- ratio.sum.stat
        pvalue.fun <- function() -1
        stat.list <- FALSE
    } else if (stat.type == 5) {
        stat.fun <- ratio.sum.stat2
        pvalue.fun <- function() -1
        stat.list <- FALSE
    } else if (stat.type == 6) {
        stat.fun <- ratio.sum.stat3
        pvalue.fun <- function() -1
        stat.list <- FALSE
    } else {
        # default
        stat.fun <- scale.stat3
        pvalue.fun <- scale.pvalue2
        stat.list <- TRUE
    }

    if (stat.list) {
        stat <- sapply(wavcoefs, stat.fun)
        pvalue <- sapply(wavcoefs, pvalue.fun)
    } else {
        stat <- stat.fun()
        pvalue <- pvalue.fun()
    }

    if (debug) {
        list(#y.up=y.up, y.lo=y.lo, 
             var.up=var.up, var.lo=var.lo, 
             Mj=Mj, #Aj=Aj, Aj.up=Aj.up, Aj.lo=Aj.lo, 
             stat=stat, pvalue=pvalue, #df.up=df.up, df.lo=df.lo, 
             reject=pvalue["d1"]<pvalue.alpha)
    } else {
        stat
    }
}


surr.wav.test <- function(y, delay, surr.N=100, q=0.333, wavelet="s8", 
                      pvalue.alpha=0.05, levels=NULL, aaft=FALSE) {
    modwt.missing.opts <- list(r=modwt.missing.r,
                               rcpp=modwt.missing.rcpp)
    modwt.missing <- modwt.missing.opts[[getOption("setarwav.mode")]]

    if (is.null(modwt.missing)) {
        stop("getOption(\"setarwav\")$modwt.missing should be either \"r\" or \"rcpp\"")
    }

    stat <- function(x) {
        wav.test(x, delay=delay, q=q, wavelet=wavelet, debug=FALSE, 
                  modwt.missing=modwt.missing,
                  levels=levels, stat.type=6) 
    }

    surr <- surrogate(y, ns=surr.N, fft=TRUE, amplitude=aaft, statistic=stat)

    y.stat <- surr$orig.statistic
    surr.data <- with(surr, statistic[!is.na(statistic)])

    # 2-sided
    #conf.int <- quantile(surr$statistic, c(pvalue.alpha/2,1-pvalue.alpha/2))
    #reject <- y.stat < conf.int[1] || y.stat > conf.int[2]

    # 1-sided
    conf.int <- quantile(surr.data, 1-pvalue.alpha, na.rm=TRUE)
    reject <- y.stat > conf.int[1]

    surr.ecdf <- ecdf(surr.data)
    pvalue <- 1 - surr.ecdf(y.stat)

    # 2-sided
    #if (pvalue > 0.5) {
    #    pvalue <- 2 * (1 - pvalue)
    #} else {
    #    pvalue <- 2 * pvalue
    #}

    # cleanup names
    names(reject) <- NULL

    list(reject=reject, stat=y.stat, conf.int=conf.int, pvalue=pvalue)

}
