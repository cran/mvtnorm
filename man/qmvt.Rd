\name{qmvt}
\alias{qmvt}
\title{ Quantiles of the Multivariate t Distribution }
\description{

Computes the equicoordinate quantile function of the multivariate t
distribution for arbitrary correlation matrices
based on an inversion of the algorithms by Genz and Bretz.

}
\usage{
qmvt(p, interval = c(-10, 10), tail = c("lower.tail", 
     "upper.tail", "both.tails"), df = 1, delta = 0, corr = NULL, 
     sigma = NULL, algorithm = GenzBretz(), ...)
}
\arguments{
  \item{p}{ probability.}
  \item{interval}{ a vector containing the end-points of the interval to be
          searched by \code{\link{uniroot}}.}
  \item{tail}{ specifies which quantiles should be computed.
               \code{lower.tail} gives the quantile \eqn{x} for which 
               \eqn{P[X \le x] = p}, \code{upper.tail} gives \eqn{x} with 
               \eqn{P[X > x] = p} and 
               \code{both.tails} leads to \eqn{x} 
               with \eqn{P[-x \le X \le x] = p}.}
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
  \item{...}{ additional parameters to be passed to                                            
              \code{\link{uniroot}}.}
}
\details{

  Only equicoordinate quantiles are computed, i.e., the quantiles in each
  dimension coincide. Currently, the distribution function is inverted by 
  using the
  \code{\link{uniroot}} function which may result in limited accuracy of the 
  quantiles. 
}
\value{
  A list with four components: \code{quantile} and \code{f.quantile} 
  give the location of the quantile and the value of the function
  evaluated at that point. \code{iter} and \code{estim.prec} give the number 
  of iterations used and an approximate estimated precision from
  \code{\link{uniroot}}.
}
\seealso{\code{\link{pmvnorm}}, \code{\link{qmvt}}}
\examples{
qmvt(0.95, df = 16, tail = "both")
}
\keyword{distribution}
