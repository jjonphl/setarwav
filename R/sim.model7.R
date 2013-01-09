sim.model7 <-
function(n, sd=1, startup=100, full=FALSE) {
    setar.sim.simple(thresholds=0, model=list(c(0, 0.6), c(0, 0.6)), 
                      n=n, delay=1, sd=c(sd, 1.5*sd), startup=startup)
}
