sim.model6 <-
function(n, sd=1, startup=500, full=FALSE) {
    setar.sim.simple(thresholds=1, model=list(c(2, 0.5), c(0.5, -0.4)), 
                     n=n, delay=1, sd=sd, startup=startup)
}
