tsaypet <- function(y, maxp, maxd) {
    nobs <- length(y)
    pvtsay <- 2            # p-value tsay
    pvpetru <- .11         # p-value petrucelli

    for (delay in 1:maxd) {
        for (p in 1:maxp) {
            df1 <- p + 1
            h <- max(1, p+1-delay)
            b <- floor(nobs/10) + p
            X <- matrix(1, b, p+1)
            yx <- matrix(0, b, 1)
            ndx <- h:(nobs-delay)
            ytd <- y[ndx]
            ix <- sort(ytd, index.return=TRUE)$ix
            yt <- as.matrix(cbind(ndx[ix], ytd[ix]))

            j <- 0; j1 <- 0
            while (j < b) {
                j1 <- j1 + 1
                pii <- yt[j1, 1]
                k <- pii + delay
                if (k > p) {
                    j <- j+1; yx[j] <- y[k];
                    for (p1 in 2:df1) {
                        X[j,p1] <- y[k-p1+1]
                    }
                    #p1 <- 2:df1
                    #X[j,p1] <- y[k-p1+1]
                }
            }

            pm <- solve(t(X) %*% X)
            beta <- pm %*% t(X) %*% yx
            s <- sum((yx - X %*% beta)^2)
            eff <- nobs - delay - h + 1 - b
            xm <- matrix(1, p+1, 1)
            xm2 <- matrix(1, eff, p+1)
            resid <- matrix(0, eff, 1)
            z <- matrix(0, eff, 1)

            for (i in (b+1):(nobs-delay-h+1)) {
                i1 <- i - b
                div <- i - p - 1
                pii <- yt[i, 1]
                k <- pii + delay
                ym <- y[k]

                #for (p1 in 2:(p+1)) {
                #    xm[p1,1] <- y[k-p1+1]
                #}
                p1 <- 2:(p+1)
                xm[p1,1] <- y[k+1-p1]

                xm2[i1,] <- xm
                dm <- as.numeric(1 + t(xm) %*% pm %*% xm)
                km <- (pm %*% xm)/dm
                a <- ym - t(xm) %*% beta
                inter <- (xm %*% t(xm))/dm
                pm <- (diag(1,p+1) - pm %*%inter) %*% pm
                beta <- beta + km %*% (ym - t(xm) %*% beta)
                resid[i1,] <- a/sqrt(dm)
                s <- s + resid[i1,]^2
                z[i1,] <- resid[i1]/sqrt(s/div)
            }

            nb <- solve(t(xm2) %*% xm2) %*% (t(xm2) %*% resid)
            eps <- resid - (xm2 %*% nb)
            sumsqres <- sum(resid^2)
            sumsqeps <- sum(eps^2)
            df2 <- nobs - delay - b - p - h
            F <- ((sumsqres - sumsqeps)*df2)/(sumsqeps*df1);
            pvtsay <- round(1 - pf(F,df1, df2), 4)

            nstar <- eff
            revcumz <- matrix(0, nstar, 1)
            revcumz[1] <- z[nstar]

            for (i in 2:nstar) {
                subs <- nstar + 1 - i
                revcumz[i,] <- revcumz[i-1,] + z[subs,]
            }
            #i <- 2:nstar
            #revcumz[i] <- revcumz[i-1] + z[nstar+1-i]

            u = 0
            alpha <- 0.001
            while (alpha <= 0.1 && u != 1) {
                i <- 1
                while (i <= nstar && u != 1) {
                    a <- abs(revcumz[i,])
                    slope <- 2*sqrt(log(2/alpha)/(2*nstar))
                    intrcpt <- 0.5*sqrt(nstar*log(2/alpha)/2)
                    a1 <- slope*i + intrcpt
                    if (a > a1) {
                        pvpetru <- alpha
                        u <- 1
                    }
                    i <- i+1
                }
                alpha <- alpha + 0.0005
            }
        }  # p
    } # d
    list(delay=maxd, p=maxp, tsayFdf=c(df1,df2), pvtsay=pvtsay, pvp=pvpetru)
} # function
