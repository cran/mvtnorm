\name{pmvt}
\alias{pmvt}
\alias{rmvt}
\title{ Multivariate t Distribution }
\description{

Computes the the distribution function of the multivariate t distribution 
for arbitrary limits, degrees of freedom and correlation matrices 
based on algorithms by Genz and Bretz.

}
\usage{
pmvt(lower=-Inf, upper=Inf, delta=rep(0, length(lower)),
     df=1, corr=NULL, sigma=NULL, maxpts = 25000, abseps = 0.001,
     releps = 0)
rmvt(n, sigma=diag(2), df=1)
}
\arguments{
  \item{lower}{ the vector of lower limits of length n.}
  \item{upper}{ the vector of upper limits of length n.}
  \item{delta}{ the vector of noncentrality parameters of length n.}
  \item{df}{ degree of freedom as integer.}
  \item{corr}{ the correlation matrix of dimension n.}
  \item{sigma}{ the covariance matrix of dimension n. Either \code{corr} or
                \code{sigma} can be specified. If \code{sigma} is given, the  
                problem is standardized. If neither \code{corr} nor
                \code{sigma} is given, the identity matrix is used 
                for \code{sigma}. }
  \item{maxpts}{ maximum number of function values as integer. }
  \item{abseps}{ absolute error tolerance as double. }
  \item{releps}{ relative error tolerance as double. }
  \item{n}{ number of observations.}
}

\details{

This program involves the computation of central and noncentral
multivariate t-probabilities with arbitrary correlation matrices.
It involves both the computation of singular and nonsingular
probabilities. The methodology is described in
Genz and Bretz (1999, 2001).

For a given correlation matrix \code{corr}, for short \eqn{A}, say, 
(which has to be  positive semi-definite) and 
degrees of freedom \code{df} the following 
values are numerically evaluated

\deqn{I = K \int s^{df-1} \exp(-s^2/2) \Phi(s \cdot lower/\sqrt{df}-delta,
s \cdot upper/\sqrt{df}-delta) ds }

where \eqn{\Phi(a,b) = K^\prime \int_a^b \exp(-x^\prime Ax/2) dx} is the
multivariate normal distribution, \eqn{K^\prime = 1/\sqrt{det(A)(2\pi)^m}}
and \eqn{K = 2^{1-df/2} / Gamma(df/2)} are constants and the (single)
integral of \eqn{I} goes from 0 to +Inf.

Note that both \code{-Inf} and \code{+Inf} may be specified in the lower and upper integral
limits. Randomized quasi-Monte Carlo methods are used for the computations.

Univariate problems are passed to \code{\link{pt}}.

Further information can be obtained from the quoted articles,
which can be downloaded (together with additional material
and additional codes) from the websites 
\url{http://www.bioinf.uni-hannover.de/~bretz/} 
and \url{http://www.sci.wsu.edu/math/faculty/genz/homepage}.

\code{rmvt} is a wrapper to \code{\link{rmvnorm}} for random number
generation.

}

\value{
The distribution function is returned with attributes
  \item{error}{estimated absolute error and}
  \item{msg}{status messages.}
}

\references{ 

Genz, A. and Bretz, F. (1999), Numerical computation of multivariate
t-probabilities with application to power calculation of multiple
contrasts. \emph{Journal of Statistical Computation and Simulation},
\bold{63}, 361--378.

Genz, A. and Bretz, F. (2001), Methods for the computation of multivariate
t-probabilities. \emph{(submitted)}

Edwards D. and Berry, Jack J. (1987), The efficiency of simulation-based 
multiple comparisons. \emph{Biometrics}, \bold{43}, 913--928

}

\author{ Fortran Code by Alan Genz <AlanGenz@wsu.edu> and Frank Bretz
<bretz@ifgb.uni-hannover.de>, \R port by Torsten Hothorn
<Torsten.Hothorn@rzmail.uni-erlangen.de> }

\examples{

n <- 5
lower <- -1
upper <- 3
df <- 4
corr <- diag(5)
corr[lower.tri(corr)] <- 0.5
delta <- rep(0, 5)
prob <- pmvt(lower=lower, upper=upper, delta=delta, df=df, corr=corr)
print(prob)

pmvt(lower=-Inf, upper=3, df = 3, sigma = 1) == pt(3, 3)

# Example from R News paper (original by Edwards and Berry, 1987)

n <- c(26, 24, 20, 33, 32)
V <- diag(1/n)
df <- 130
C <- c(1,1,1,0,0,-1,0,0,1,0,0,-1,0,0,1,0,0,0,-1,-1,0,0,-1,0,0)
C <- matrix(C, ncol=5)
cv <- C \%*\% V \%*\% t(C)
cr <- matrix(rep(0, ncol(cv)^2), ncol=ncol(cv))
for (i in 1:5) {
  for (j in 1:5) {
    cr[i,j] <- cv[i,j]/sqrt(cv[i,i]*cv[j,j] )
  }
}
delta <- rep(0,5)

myfct <- function(q, alpha) {
  lower <- rep(-q, ncol(cv))
  upper <- rep(q, ncol(cv))
  pmvt(lower=lower, upper=upper, delta=delta, df=df, corr=cr, abseps=0.0001) - alpha
}

round(uniroot(myfct, lower=1, upper=5, alpha=0.95)$root, 3)

# compare pmvt and pmvnorm for large df:

a <- pmvnorm(lower=-Inf, upper=1, mean=rep(0, 5), corr=diag(5))
b <- pmvt(lower=-Inf, upper=1, delta=rep(0, 5), df=rep(300,5),
          corr=diag(5))
a
b

stopifnot(round(a, 2) == round(b, 2))

# correlation and covariance matrix

a <- pmvt(lower=-Inf, upper=2, delta=rep(0,5), df=3,
          sigma = diag(5)*2)
b <- pmvt(lower=-Inf, upper=2/sqrt(2), delta=rep(0,5),
          df=3, corr=diag(5))  
attributes(a) <- NULL
attributes(b) <- NULL
a
b
stopifnot(all.equal(round(a,3) , round(b, 3)))

a <- pmvt(0, 1,df=10)
attributes(a) <- NULL
b <- pt(1, df=10) - pt(0, df=10)
stopifnot(all.equal(round(a,10) , round(b, 10)))

rmvt(10, sigma=diag(10))

}

\keyword{distribution}
