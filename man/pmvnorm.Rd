\name{pmvnorm}
\alias{pmvnorm}
\title{ Multivariate Normal Distribution }
\description{

Computes the distribution function of the multivariate normal 
distribution for arbitrary limits and correlation matrices 
based on algorithms by Genz and Bretz.

}
\usage{
pmvnorm(lower, upper, mean, corr, maxpts = 25000, abseps =0.001, releps = 0)
}
\arguments{
  \item{lower}{ the vector of lower limits of length n.}
  \item{upper}{ the vector of upper limits of length n.}
  \item{mean}{ the mean vector of length n.}
  \item{corr}{ the correlation matrix of dimension n.}
  \item{maxpts}{ maximum number of function values as integer. }
  \item{abseps}{ absolute error tolerance as double. }
  \item{releps}{ relative error tolerance as double. }
}
\details{

This program involves the computation of 
multivariate normal-probabilities with arbitrary correlation matrices.
It involves both the computation of singular and nonsingular 
probabilities. The methodology is described in
Genz (1992, 1993).

Note that both \code{-Inf} and \code{+Inf} may be specified in \code{lower} and
\code{upper}. For more details see \code{\link{pmvt}}. 

The mvn case is treated as a special case of \code{\link{pmvt}} with \code{df=0}. 

Univariate problems are passed to \code{\link{pnorm}}.

Multivariate normal density and random numbers are available using
\code{\link[e1071]{dmvnorm}} and \code{\link[MASS]{mvrnorm}}.  
 
}
\value{
A list with the following components:
  \item{value}{estimated integral value.}
  \item{error}{estimated absolute error.}
  \item{msg}{status messages.}
}
\references{ 

Genz, A. (1992). Numerical computation of multivariate normal probabilities.
\emph{Journal of Computational and Graphical Statistics}, \bold{1}, 141--150 

Genz, A. (1993). Comparison of methods for the computation of multivariate
normal probabilities. \emph{Computing Science and Statistics}, \bold{25},
400--405

}

\author{ Fortran Code by Alan Genz <AlanGenz@wsu.edu> and Frank Bretz
<bretz@ifgb.uni-hannover.de>, \R port by Torsten Hothorn
<Torsten.Hothorn@rzmail.uni-erlangen.de> }

\examples{

n <- 5
mean <- rep(0, 5)
lower <- rep(-1, 5)
upper <- rep(3, 5)
corr <- diag(5)
corr[lower.tri(corr)] <- 0.5
corr[upper.tri(corr)] <- 0.5
prob <- pmvnorm(lower, upper, mean, corr)
print(prob)

stopifnot(pmvnorm(lower=-Inf, upper=3, mean=0, corr=1)$value == pnorm(3))

a <- pmvnorm(lower=rep(-Inf,2),upper=c(.3,.5),mean=c(2,4),diag(2))$value

stopifnot(round(a,16) == round(prod(pnorm(c(.3,.5),c(2,4))),16))

# Example from R News paper (original by Genz, 1992):

m <- 3
sigma <- diag(3)
sigma[2,1] <- 3/5
sigma[3,1] <- 1/3
sigma[3,2] <- 11/15
pmvnorm(lower=rep(-Inf, m), upper=c(1,4,2), mean=rep(0, m), sigma)

}

\keyword{distribution}
