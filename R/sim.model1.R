sim.model1 <-
function(n) {
    arima.sim(n, model=list(ar=c(0.4, -0.3)))
}
