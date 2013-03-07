\name{sim.model*}
\alias{sim.model1}
\alias{sim.model2}
\alias{sim.model3}
\alias{sim.model4}
\alias{sim.model5}
\alias{sim.model6}
\alias{sim.model7}
\alias{sim.tsay}
\title{
Simulate a realization of a specific time series models of length n.
}
\description{
These methods are used to generate realizations of a set of time series 
models for use in evaluating the power of the statistical tests in this
package.
}
\usage{
sim.model1(n)
sim.model2(n)
sim.model3(n)
sim.model4(n)
sim.model5(n)
sim.model6(n)
sim.model7(n)
sim.tsay(n)
}
\arguments{
  \item{n}{Length of returned realization.}
}
\details{
The default burn-in period is 100, except for \code{sim.model1} and 
\code{sim.model2} uses burn which uses the burn-in period computed
by \code{arima.sim}
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
function (n) 
{
    arima.sim(n, model = list(ar = c(0.4, -0.3)))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
