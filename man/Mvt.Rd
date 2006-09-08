\name{Mvt}
\alias{dmvt}
\alias{rmvt}
\title{The Multivariate t Distribution}
\description{
  These functions provide information about the multivariate t
  distribution with non-centrality parameter \code{delta}, covariance matrix
  \code{sigma} and degrees of freedom \code{df}. 
  \code{dmvt} gives the density and \code{rmvt}
  generates random deviates.
}
\usage{
rmvt(n, sigma = diag(2), df = 1)
dmvt(x, delta, sigma, df = 1, log = TRUE)
}
\arguments{
 \item{x}{Vector or matrix of quantiles. If \code{x} is a matrix, each
   row is taken to be a quantile.}
 \item{n}{Number of observations.}
 \item{delta}{ the vector of noncentrality parameters of length n.}
 \item{sigma}{Covariance matrix, default is \code{diag(ncol(x))}.}
 \item{df}{ degree of freedom as integer.}
 \item{log}{Logical; if \code{TRUE}, densities d are given as log(d).}
}
\seealso{\code{\link{pmvt}} and \code{\link{qmvt}}}
\examples{

  dmvt(x=c(0,0), sigma = diag(2))
  x <- rmvt(n=100, sigma = diag(2), df = 3)
  plot(x)

}
\keyword{distribution}
\keyword{multivariate}
