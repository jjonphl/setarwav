sim.tsay <-
function(n, phi=-2, sd=1, startup=500) {
    setar.sim.simple(thresholds=0, model=list(c(0, 0.5), c(0, phi)), n=n,
                     delay=1, startup=startup)
}
