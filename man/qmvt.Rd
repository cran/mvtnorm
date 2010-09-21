\name{qmvt}
\alias{qmvt}
\title{ Quantiles of the Multivariate t Distribution }
\description{

Computes the equicoordinate quantile function of the multivariate t
distribution for arbitrary correlation matrices
based on inversion of qmvt.

}
\usage{
qmvt(p, interval = NULL, tail = c("lower.tail", 
     "upper.tail", "both.tails"), df = 1, delta = 0, corr = NULL, 
     sigma = NULL, algorithm = GenzBretz(), 
     type = c("Kshirsagar", "shifted"), ...)
}
\arguments{
  \item{p}{ probability.}
  \item{interval}{ optional, a vector containing the end-points of the interval to be
          searched by \code{\link{uniroot}}.}
  \item{tail}{ specifies which quantiles should be computed.
               \code{lower.tail} gives the quantile \eqn{x} for which 
               \eqn{P[X \le x] = p}, \code{upper.tail} gives \eqn{x} with 
               \eqn{P[X > x] = p} and 
               \code{both.tails} leads to \eqn{x} 
               with \eqn{P[-x \le X \le x] = p}.}
  \item{delta}{ the vector of noncentrality parameters of length n, for
   \code{type = "shifted"} delta specifies the mode.}
  \item{df}{ degree of freedom as integer. Normal quantiles are computed for \code{df = 0}.}
  \item{corr}{ the correlation matrix of dimension n.}
  \item{sigma}{ the covariance matrix of dimension n. Either \code{corr} or
                \code{sigma} can be specified. If \code{sigma} is given, the
                problem is standardized. If neither \code{corr} nor
                \code{sigma} is given, the identity matrix is used
                for \code{sigma}. }
  \item{algorithm}{ an object of class \code{\link{GenzBretz}} or
                    \code{\link{TVPACK}} defining the
                    hyper parameters of this algorithm.}
  \item{type}{ type of the noncentral multivariate t distribution
              to be computed. \code{type = "Kshirsagar"} corresponds
              to formula (1.4) in Genz and Bretz (2009) (see also
	      Chapter 5.1 in Kotz and Nadarajah (2004)) and
	      \code{type = "shifted"} corresponds to the formula before
	      formula (1.4) in Genz and Bretz (2009)
	      (see also formula (1.1) in Kotz and Nadarajah (2004)). }
  \item{...}{ additional parameters to be passed to                                            
              \code{\link{GenzBretz}}.}
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
\seealso{\code{\link{pmvnorm}}, \code{\link{qmvnorm}}}
\examples{
qmvt(0.95, df = 16, tail = "both")
}
\keyword{distribution}
