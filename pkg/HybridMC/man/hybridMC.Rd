\name{hybridMC}
\Rdversion{1.1}
\alias{hybridMC}
\title{
Sample from a (log)density using hybrid Monte Carlo
}
\description{
\code{hybridMC()} samples from a (log) joint density function defined up to a proportionality constant. It uses the 
multipoint hybrid Monte Carlo methods described by Liu (2001). The function supports a full range of tweaking options to
minimize autocorrlation.
}
\usage{
hybridMC(y.start, n.samp=1, logDens, dLogDens, epsilon, LFsteps, compWeights=NULL, MPwidth=1, MPweights=NULL, progress=0 ,...)
}
\arguments{
  \item{y.start}{
A vector of starting values
}
  \item{n.samp}{
The number of samples to draw. Each previous value is used as the starting value for the next sample
}
  \item{logDens}{
The log-density function (up to an additive constant) from which you would like to sample. \code{logDens} must return a single value
}
  \item{dLogDens}{
The function giving the derivative of the log-density function with respect to each variable. \code{dLogDens} must return a vector of the same length as \code{y.start}
}
  \item{epsilon}{
Either a single positive value giving the size of the time simulation steps, or a vector of length 2 giving the lower and upper bounds of the interval from which \code{epsilon} is uniformly sampled
}
  \item{LFsteps}{
An integer giving the number of leapfrog simulation steps to do
}
  \item{compWeights}{
	The ``masses'' of the dimensions. Must be either a single numeric value or a vector of the same length as \code{y.start}
}
  \item{MPwidth}{
The (integer) size of the multipoint window. The default is 1, meaning no multipoint
\code{MPwidth} must be less than or equal to \code{LFsteps}.
}
  \item{MPweights}{
A vector of length \code{MPwidth} of constants, used to weight the proposal values within the multipoint window. If only a single value is passed, the values are weighted evenly
}
  \item{progress}{
An integer giving the number of samples between updates to a text progress bar. If \code{progress=0}, no progress bar is displayed 
}

\item{\dots}{
     Arguments to pass to \code{logDens} and \code{dLogDens}
}
}
\details{

The density should have support of the whole real line for every dimension. The quality of the samples can be improved by tweaking \code{epsilon}, \code{LFsteps}, \code{MPwidth}, and to to a lesser extent, \code{compWeights} and \code{MPweights}. 

If you use the progress bar, keep in mind that updating the progress bar takes a little time. Too many updates will slow down your sampling, so pick a reasonable value for \code{progress}.
}

\value{
\code{hybridMC()} returns an object of type \code{mcmc} (from the \code{coda} package). Each row is a sample from the joint density.
}

\references{
Liu (2001) "Monte Carlo strategies in scientific computing"
}
\author{
Richard D. Morey
}
\examples{

#### Example 1: Jointly sample from two independent double exponentials

## log density function for double exponential
L = function(x)
-sum(abs(x))

dL = function(x)
-sign(x)

startVal=c(-1,1)
samples = hybridMC(y.start=startVal, n.samp=5000, logDens=L, dLogDens=dL, epsilon=.2, LFsteps=10)

# Plot the MCMC chains and densities
plot(samples)

# Plot a histogram of the first variable, with true density
hist(samples[,1],freq=FALSE,breaks=50)
x = seq(-5,5,len=100)
lines(x,.5*dexp(abs(x)),col="red") #true density in red

# Autocorrelation function
acf(samples)

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{distribution}% at least one, from doc/KEYWORDS
\keyword{multivariate}% __ONLY ONE__ keyword per line
\keyword{misc}
