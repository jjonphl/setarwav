
spect.sas <- function(y, maxd, q=0.333) {
    ss <- quantile(y, prob=c(q, 1-q), type=4)
    q1 <- ss[1]; q3 <- ss[2]
    nobs <- length(y)

    Mset   <- seq.int(10, 3, by=-1)
    Mset.l <- length(Mset)
    nr     <- Mset.l * maxd
    b <- data.frame(M=integer(0) , delay=integer(0), chisq=numeric(0),
                    df=numeric(0), pvalue=numeric(0))
    f <- data.frame(M=integer(0), delay=integer(0), chisq=numeric(0),
                    df=numeric(0), pvalue=numeric(0))
    p <- data.frame(M=integer(0), delay=integer(0), chisq=numeric(0),
                    df=numeric(0), pvalue=numeric(0))
    for (del in 1:maxd) {
        eff <- nobs - del
        startx <- del + 1
        effnob <- nobs - del
        width <- floor(sqrt(effnob))
        v <- y[startx:nobs]

        upper <- as.numeric(y[1:eff] > q1)    # indicator if y > q1
        lower <- as.numeric(y[1:eff] < q3)    # indicator if y < q3

        zl <- v * lower; zh <- v * upper;
        zlsum <- sum(zl);  ilsum <- sum(lower); zlmean <- zlsum/ilsum;
        zhsum <- sum(zh);  ihsum <- sum(upper); zhmean <- zhsum/ihsum;

        zl <- zl - zlmean; zh <- zh - zhmean;   # center data, i.e. mean 0
        zl <- zl * lower; zh <- zh * upper;     # reset to 0

        # compute autocovariances up to lag 26
        zlaut <- acf(zl, 27-1, type="cov", demean=F, plot=F)$acf
        zhaut <- acf(zh, 27-1, type="cov", demean=F, plot=F)$acf
        ilaut <- acf(lower, 27-1, type="cov", demean=F, plot=F)$acf
        ihaut <- acf(upper, 27-1, type="cov", demean=F, plot=F)$acf

        startM <- max(width-2, 3)     # not used??
        endM   <- startM + 5

        pchi <- 0
        #bspect <- list()

        for (M in Mset) {
            em <- M
            M1 <- M + 1

            #-----------------------------------------------------
            # PART 1 - bartlett-priestly
            #-----------------------------------------------------
            endj <- floor(em/.775)
            delw <- (.775*pi / em)
            #bspect[[M]] <- rep(0, endj)    # JON: save computed spect

            # bartlett-priestly
            r <- 1:M
            sm <- r * pi/M
            wt <- 3*((sin(sm)/sm) - cos(sm))/(sm^2)

            for (d in 1) {
                df1 <- 0
                chi1 <- 0
                pchi <- 0
                lvar <- 1/(effnob * ilaut[1])
                hvar <- 1/(effnob * ihaut[1])

                for (j in 1:endj) {  # iterate through freq
                    w <- delw * j
                    lspectrum <- zlaut[1] / ilaut[1]
                    hspectrum <- zhaut[1] / ihaut[1]

                    for (k2 in 2:M1) {  # iterate through time
                        k1 <- k2 - 1
                        if (ilaut[k2] != 0) {
                            lspectrum <- lspectrum +
                               2*cos(w*k1)*wt[k1]*zlaut[k2]/ilaut[k2]
                        }
                        if (ihaut[k2] != 0) {
                            hspectrum <- hspectrum +
                               2*cos(w*k1)*wt[k1]*zhaut[k2]/ihaut[k2]
                        }
                    }

                    #bspect[[M]][j] <- lspectrum # JON: save computed spect

                    if (lspectrum * hspectrum > 0) {
                        xn <- log(lspectrum / hspectrum)
                        chi1 <- chi1 + xn^2
                        df1 <- df1 + 1
                    }
                }  # j

                for (j2 in 2:M1) {
                    il <- ilaut[j2]
                    ih <- ihaut[j2]
                    j1 <- j2-1
                    if (il != 0) lvar <- lvar + 2*(wt[j1]^2) / (effnob*il)
                    if (ih != 0) hvar <- hvar + 2*(wt[j1]^2) / (effnob*ih)
                }

                if (isTRUE(all.equal(endj, em/0.775))) {
                    df1 <- df1 - 1
                    chi1 <- chi1 - xn^2
                }

                chi1 <- chi1 / (lvar+hvar)

                if (df1 > 0) {
                    pchi <- round(1 - pchisq(chi1, df1), 7)
                }
            } # d

            b <- rbind(b, data.frame(M=em,delay=del,chisq=chi1,df=df1, pchi=pchi))

            #-----------------------------------------------------
            # PART 2 - fejer
            #-----------------------------------------------------

            endj <- floor(pi*em/2.45)
            delw <- (2.45/em)

            r <- 1:em
            wt <- 1 - r/em    # window

            for (d in 1) {
                df1 <- 0;
                chi1 <- 0
                pchi <- 0
                lvar <- 1/(effnob * ilaut[1])
                hvar <- 1/(effnob * ihaut[1])

                for (j in 1:endj) {
                    w <- delw * j
                    lspectrum <- zlaut[1] / ilaut[1]
                    hspectrum <- zhaut[1] / ihaut[1]

                    for (k2 in 2:M1) {
                        k1 <- k2 - 1
                        if (ilaut[k2] != 0) {
                            lspectrum <- lspectrum + 
                                2*cos(w*k1)*wt[k1]*zlaut[k2]/ilaut[k2]
                        }
                        if (ihaut[k2] != 0) {
                            hspectrum <- hspectrum + 
                                2*cos(w*k1)*wt[k1]*zhaut[k2]/ihaut[k2]
                        }
                    } # k2

                    if (lspectrum * hspectrum > 0) {
                        xn <- log(lspectrum / hspectrum) 
                        chi1 <- chi1 + xn^2
                        df1 <- df1 + 1
                    }
                } # j

                for (j2 in 2:M1) {
                    il <- ilaut[j2]
                    ih <- ihaut[j2]
                    j1 <- j2 - 1
                    if (il != 0) lvar <- lvar + 2*(wt[j1]^2) / (effnob*il)
                    if (ih != 0) hvar <- hvar + 2*(wt[j1]^2) / (effnob*ih)
                }

                if (isTRUE(all.equal(endj, pi*em/2.45))) {
                    df1 <- df1 - 1
                    chi1 <- chi1 - xn^2
                }
                chi1 <- chi1 / (lvar+hvar)

                if (df1 > 0) {
                    pchi <- round(1 - pchisq(chi1, df1), 7)
                } else {
                    pchi <- 0
                }
            }  # d

            f <- rbind(f, data.frame(M=em,delay=del,chisq=chi1,df=df1,pchi=pchi))

            #-----------------------------------------------------
            # PART 3 - parzen
            #-----------------------------------------------------
            em <- round(em*2.5)
            M1 <- em + 1
            endj <- floor(pi*em/6)
            delw <- 6/em

            r <- 1:em
            sm <- r/em
            wt <- ifelse(sm < 0.5, 1 - 6*sm^2 + 6*sm^3, 2*(1-sm)^3)

            for (d in 1) {
                df <- 0
                ch1 <- 0
                pchi <- 0
                lvar <- 1/(effnob*ilaut[1])
                hvar <- 1/(effnob*ihaut[1])

                for (j in 1:endj) {
                    w <- delw*j
                    lspectrum <- zlaut[1] / ilaut[1]
                    hspectrum <- zhaut[1] / ihaut[1]

                    for (k2 in 2:M1) {
                        k1 <- k2 - 1
                        if (ilaut[k2] != 0) {
                            lspectrum <- lspectrum +
                                2*cos(w*k1)*wt[k1]*zlaut[k2]/ilaut[k2]
                        }
                        if (ihaut[k2] != 0) {
                            hspectrum <- hspectrum +
                                2*cos(w*k1)*wt[k1]*zhaut[k2]/ihaut[k2]
                        }
                    }

                    if (lspectrum * hspectrum > 0) {
                        xn <- log(lspectrum/hspectrum)
                        chi1 <- chi1 + xn**2
                        df1 <- df1 + 1
                    }
                } # j

                for (j2 in 2:M1) {
                    il <- ilaut[j2]
                    ih <- ihaut[j2]
                    j1 <- j2 - 1
                    if (il != 0) lvar <- lvar+2*(wt[j1]^2)/(effnob*il)
                    if (ih != 0) hvar <- hvar+2*(wt[j1]^2)/(effnob*ih)
                }

                if (TRUE) {
                    df1 <- df1 - 1
                    chi1 <- chi1 - xn^2
                }
                chi1 <- chi1 / (lvar+hvar)

                if (df1 > 0) {
                    pchi <- round(1 - pchisq(chi1, df1), 7)
                } else {
                    pchi <- 0
                }
            } # d

            p <- rbind(p, data.frame(M=em,delay=del,chisq=chi1,df=df1,pchi=pchi))
        } # M
    } # del

    list(bartlett=b, fejer=f, parzen=p)
}


