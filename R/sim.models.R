# linear AR(2): X[n] = 0.4*X[n-1] - 0.3X[n-2] + e[n]
sim.model01 <- function(n) {
    arima.sim(n, model=list(ar=c(0.4, -0.3)))
}

# linear MA(2) - X[n] = e[n] - 0.4e[n-1] - 0.3e[n-2]
sim.model02 <- function(n) {
    arima.sim(n, model=list(ar=c(-0.4, 0.3)))
}

# EXPAR(1) - X[n] = (0.3 - 0.8exp(-X[n-1]^2))*X[n-1] + e[n]
sim.model03 <- function(n, startup=100) {
    expar.sim(n, model=list(a=0.3,b=-0.8,c=0,d=1), d=1, sd=1, startup=100)
}

# Bilinear - X[n] = 0.5 - 0.4X[n-1] + 0.4*X[n-1]*e[n-1] + e[n]
sim.model04 <- function(n, startup=100) {
    bilinear.sim(n, model=list(a0=0.5,ar=-0.4,bl=0.4), 
                 sd=1, startup=startup)
}

# Bilinear - X[n] = 0.4X[n-1] - 0.3X[n-2] + 0.5X[n-1]*e[n-1] + 0.8e[n-1] + e[n]
sim.model05 <- function(n, startup=100) {
    bilinear.sim(n, model=list(ar=c(0.4,-0.3),ma=0.8,bl=0.5),
                 sd=1, startup=startup)
}

# SETAR - X[n] = | X[n-1] < 1: 2 + 0.5*X[n-1] + e[n]
#                | else: 0.5 - 0.4*X[n-1] + e[n]
sim.model06 <- function(n, sd=1, startup=100, full=FALSE) {
    setar.sim.simple(thresholds=1, model=list(c(2, 0.5), c(0.5, -0.4)), 
                     n=n, delay=1, sd=sd, startup=startup)
}

# SETAR  - X[n] = | X[n-1] < 0: 0.6*X[n-1] + e[n]
#                 | else: 0.6*X[n-1] + 1.5*e[n]
sim.model07 <- function(n, sd=1, startup=100, full=FALSE) {
    setar.sim.simple(thresholds=0, model=list(c(0, 0.6), c(0, 0.6)), 
                      n=n, delay=1, sd=c(sd, 1.5*sd), startup=startup)
}

# SETAR  - X[n] = | X[n-2] < 0: 1.22*X[n-1] - 0.44*X[n-2] + e[n]
#                 | else: -0.25*X[n-1] - 0.05*X[n-2] + e[n]
sim.model08 <- function(n, sd=1, startup=100) {
    setar.sim.simple(thresholds=0, 
                     model=list(c(0, 1.22, -0.44), c(0, -0.25, -0.05)), 
                     n=n, delay=2, sd=1, startup=startup) 
}

# SETAR  - X[n] = | X[n-2] < 0: 1.22*X[n-1] - 0.44*X[n-2] + e[n]
#                 | else: 1.10*X[n-1] - 0.5*X[n-2] + e[n]
sim.model09 <- function(n, sd=1, startup=100) {
    setar.sim.simple(thresholds=0, 
                     model=list(c(0, 1.22, -0.44), c(0, 1.10, -0.5)), 
                     n=n, delay=2, sd=1, startup=startup)
}

# SETAR - X[n] = | X[n-3] < -1: -0.03*X[n-1] + 0.73*X[n-2] - 0.09*X[n-3] + e[n]
#                | X[n-3] \in [-1,1]: 1.1*X[n-1] - 0.57*X[n-2] + e[n]
#                | else: 0.52*X[n-1] - 0.45*X[n-2] + 0.14*X[n-3] + e[n]
sim.model10 <- function(n, sd=1, startup=100) {
    setar.sim.simple(thresholds=c(-1,1), 
                     model=list(c(0, -0.03, 0.73, -0.09), 
                                c(0, 1.1, -0.57),
                                c(0, 0.52, -0.45, 0.14)), 
                     n=n, delay=3, sd=1, startup=startup)
}

# SETAR - X[n] = | X[n-1] > 0: phi*X[n-1] + e[n]
#                | else: 0.5*X[n-1] + e[n]
# phi \in {-2, -1, -0.5, 0, 0.5}
sim.tsay <- function(n, phi=-2, sd=1, startup=500) {
    setar.sim.simple(thresholds=0, model=list(c(0, 0.5), c(0, phi)), n=n,
                     delay=1, startup=startup)
}


