# linear AR(2): X[n] = 0.4*X[n-1] - 0.3X[n-2] + e[n]
sim.model1 <- function(n) {
    arima.sim(n, model=list(ar=c(0.4, -0.3)))
}

# linear MA(2) - X[n] = e[n] - 0.4e[n-1] - 0.3e[n-2]
sim.model2 <- function(n) {
    arima.sim(n, model=list(ar=c(-0.4, 0.3)))
}

# EXPAR(1) - X[n] = (0.3 - 0.8exp(-X[n-1]^2))*X[n-1] + e[n]
sim.model3 <- function(n) {
}

# Bilinear - X[n] = 0.5 - 0.4X[n-1] + 0.4*X[n-1]*e[n-1] + e[n]
sim.model4 <- function(n) {
}

# Bilinear - X[n] = 0.4X[n-1] - 0.3X[n-2] + 0.5X[n-1]*e[n-1] + 0.8e[n-1] + e[n]
sim.model5 <- function(n) {
}

# SETAR - X[n] = | X[n-1] < 1: 2 + 0.5*X[n-1]
#                | else: 0.5 - 0.4*X[n-1] + e[n]
sim.model6 <- function(n, sd=1, startup=500, full=FALSE) {
    setar.sim.simple(thresholds=1, model=list(c(2, 0.5), c(0.5, -0.4)), 
                     n=n, delay=1, sd=sd, startup=startup)
}

# SETAR  - X[n] = | X[n-1] > 0: 0.6*X[n-1] + e[n]
#                 | else: 0.6*X[n-1] + 1.5*e[n]
sim.model7 <- function(n, sd=1, startup=100, full=FALSE) {
    setar.sim.simple(thresholds=0, model=list(c(0, 0.6), c(0, 0.6)), 
                      n=n, delay=1, sd=c(sd, 1.5*sd), startup=startup)
}

# phi \in {-2, -1, -0.5, 0, 0.5}
sim.tsay <- function(n, phi=-2, sd=1, startup=500) {
    setar.sim.simple(thresholds=0, model=list(c(0, 0.5), c(0, phi)), n=n,
                     delay=1, startup=startup)
}


