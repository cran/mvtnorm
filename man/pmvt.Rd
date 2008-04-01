\name{pmvt}
\alias{pmvt}
\title{ Multivariate t Distribution }
\description{

  Computes the the distribution function of the multivariate t distribution 
  for arbitrary limits, degrees of freedom and correlation matrices 
  based on algorithms by Genz and Bretz.

}
\usage{
pmvt(lower=-Inf, upper=Inf, delta=rep(0, length(lower)),
     df=1, corr=NULL, sigma=NULL, algorithm = GenzBretz(), ...)
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
  \item{algorithm}{ an object of class \code{\link{GenzBretz}} defining the
                    hyper parameters of this algorithm.}
  \item{...}{}
}

\details{

This program involves the computation of central and noncentral
multivariate t-probabilities with arbitrary correlation matrices.
It involves both the computation of singular and nonsingular
probabilities. The methodology is described in
Genz and Bretz (1999, 2002).

For a given correlation matrix \code{corr}, for short \eqn{A}, say, 
(which has to be  positive semi-definite) and 
degrees of freedom \eqn{\nu} the following 
values are numerically evaluated

\deqn{I = 2^{1-\nu/2} / \Gamma(\nu/2) \int_0^\infty s^{\nu-1} \exp(-s^2/2) \Phi(s \cdot lower/\sqrt{\nu} - \delta,
s \cdot upper/\sqrt{\nu} - \delta) \, ds }

where 

\deqn{\Phi(a,b) = (det(A)(2\pi)^m)^{-1/2} \int_a^b \exp(-x^\prime Ax/2) \, dx} 

is the multivariate normal distribution and $m$ is the number of rows of \eqn{A}.

Note that both \code{-Inf} and \code{+Inf} may be specified in 
the lower and upper integral limits in order to compute one-sided
probabilities. Randomized quasi-Monte Carlo methods are used for the 
computations.

Univariate problems are passed to \code{\link{pt}}. 
If \code{df = 0}, normal probabilities are returned. 

}

\value{
The evaluated distribution function is returned with attributes
  \item{error}{estimated absolute error and}
  \item{msg}{status messages.}
}

\references{ 

Genz, A. and Bretz, F. (1999), Numerical computation of multivariate
t-probabilities with application to power calculation of multiple
contrasts. \emph{Journal of Statistical Computation and Simulation},
\bold{63}, 361--378.

Genz, A. and Bretz, F. (2002), Methods for the computation of multivariate
t-probabilities. \emph{Journal of Computational and Graphical Statistics},
\bold{11}, 950--971.

Edwards D. and Berry, Jack J. (1987), The efficiency of simulation-based 
multiple comparisons. \emph{Biometrics}, \bold{43}, 913--928.

}

\source{
  \url{http://www.sci.wsu.edu/math/faculty/genz/homepage}
}

\seealso{\code{\link{qmvt}}}

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
### covariance matrix
cv <- C \%*\% V \%*\% t(C)
### correlation matrix
dv <- t(1/sqrt(diag(cv)))
cr <- cv * (t(dv) \%*\% dv)
delta <- rep(0,5)

myfct <- function(q, alpha) {
  lower <- rep(-q, ncol(cv))
  upper <- rep(q, ncol(cv))
  pmvt(lower=lower, upper=upper, delta=delta, df=df, 
       corr=cr, abseps=0.0001) - alpha
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

}

\keyword{distribution}
