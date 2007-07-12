\name{Mvnorm}
\alias{dmvnorm}
\alias{rmvnorm}
\title{Multivariate Normal Density and Random Numbers}
\description{
  These functions provide the density function and a random number
  generator for the multivariate normal
  distribution with mean equal to \code{mean} and covariance matrix
  \code{sigma}. 
}
\usage{
dmvnorm(x, mean, sigma, log=FALSE)
rmvnorm(n, mean, sigma, method=c("svd", "chol"))
}
\arguments{
 \item{x}{Vector or matrix of quantiles. If \code{x} is a matrix, each
   row is taken to be a quantile.}
 \item{n}{Number of observations.}
 \item{mean}{Mean vector, default is \code{rep(0, length = ncol(x))}.}
 \item{sigma}{Covariance matrix, default is \code{diag(ncol(x))}.}
 \item{log}{Logical; if \code{TRUE}, densities d are given as log(d).}
 \item{method}{Matrix decomposition used to determine the matrix root of
  \code{sigma}, possible methods are singular value decomposition
	(\code{"svd"}, default) and Cholesky decomposition (\code{"chol"}).}
} 
\details{
  SVD will work even for numerically negative definite \code{sigma}
  (with a warning), while the Cholesky decomposition is possible only
  for non-negative 
  definite \code{sigma} and issues a warning for positive
  semi-definiteness.
  Using \code{method = "chol"} may be faster and more stable for higher
  dimensional \code{sigma}.
}

\author{Friedrich Leisch <Friedrich.Leisch@R-project.org>, Fabian Scheipl}
\seealso{\code{\link{pmvnorm}}, \code{\link{rnorm}}, \code{\link{qmvnorm}}}
\examples{
dmvnorm(x=c(0,0))
dmvnorm(x=c(0,0), mean=c(1,1))

sigma <- matrix(c(4,2,2,3), ncol=2)
x <- rmvnorm(n=500, mean=c(1,2), sigma=sigma)
colMeans(x)
var(x)

x <- rmvnorm(n=500, mean=c(1,2), sigma=sigma, method="chol")
colMeans(x)
var(x)

plot(x)
}
\keyword{distribution}
\keyword{multivariate}
