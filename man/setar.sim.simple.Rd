\name{setar.sim.simple}
\alias{setar.sim.simple}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
setar.sim.simple(thresholds, model, n, delay, sd = 1, startup = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{thresholds}{
%%     ~~Describe \code{thresholds} here~~
}
  \item{model}{
%%     ~~Describe \code{model} here~~
}
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{delay}{
%%     ~~Describe \code{delay} here~~
}
  \item{sd}{
%%     ~~Describe \code{sd} here~~
}
  \item{startup}{
%%     ~~Describe \code{startup} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (thresholds, model, n, delay, sd = 1, startup = 100) 
{
    if (!is.list(model)) {
        stop("model should be a list(ar.regime1, ar.regime2, ...)")
    }
    if (length(thresholds) + 1 != length(model)) {
        stop("Length of thresholds & model params does not match.")
    }
    if (any(diff(thresholds) < 0)) {
        stop("Thresholds must be specified in ascending order")
    }
    len <- startup + n
    rs <- rep(0, len)
    thresholds <- c(-Inf, sort(thresholds), Inf)
    n.regimes <- length(thresholds) + 1
    if (length(sd) != n.regimes) 
        sd <- rep(sd, length = n.regimes)
    y <- numeric(len)
    e.list <- list()
    for (i in 1:n.regimes) {
        e.list[[i]] <- rnorm(len, sd = sd[i])
    }
    for (i in 2:len) {
        d <- max(i - delay, 1)
        regime <- which(y[d] < thresholds)[1] - 1
        coefs <- model[[regime]]
        p <- length(coefs) - 1
        idx <- seq(i - 1, max(i - p, 1))
        lags <- c(1, y[idx], rep(0, p - length(idx)))
        y[i] <- sum(lags * coefs) + e.list[[regime]][i]
    }
    return(y[(startup + 1):len])
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
