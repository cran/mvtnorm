\name{Mvnorm}
\alias{dmvnorm}
\alias{rmvnorm}
\title{The Multivariate Normal Distribution}
\description{
  These functions provide information about the multivariate normal
  distribution with mean equal to \code{mean} and covariance matrix
  \code{sigma}. \code{dmvnorm} gives the density and \code{rmvnorm}
  generates random deviates.
}
\usage{
dmvnorm(x, mean, sigma, log=FALSE)
rmvnorm(n, mean, sigma)
}
\arguments{
 \item{x}{Vector or matrix of quantiles. If \code{x} is a matrix, each
   row is taken to be a quantile.}
 \item{n}{Number of observations.}
 \item{mean}{Mean vector, default is \code{rep(0, length = ncol(x))}.}
 \item{sigma}{Covariance matrix, default is \code{diag(ncol(x))}.}
 \item{log}{Logical; if \code{TRUE}, densities d are given as log(d).}
}
\author{Friedrich Leisch <Friedrich.Leisch@ci.tuwien.ac.at> }
\seealso{\link{pmvnorm},\link{rnorm}}
\examples{
dmvnorm(x=c(0,0))
dmvnorm(x=c(0,0), mean=c(1,1))
x <- rmvnorm(n=100, mean=c(1,1))
plot(x)
}
\keyword{distribution}
\keyword{multivariate}
