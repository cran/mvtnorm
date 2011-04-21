\name{Mvt}
\alias{dmvt}
\alias{rmvt}
\title{The Multivariate t Distribution}
\description{
  These functions provide information about the multivariate t
  distribution with non-centrality parameter (or mode) \code{delta},
  covariance matrix \code{sigma} and degrees of freedom \code{df}. 
  \code{dmvt} gives the density and \code{rmvt}
  generates random deviates.
}
\usage{
rmvt(n, sigma = diag(2), df = 1, delta = rep(0, nrow(sigma)),
     type = c("shifted", "Kshirsagar"), ...)
dmvt(x, delta, sigma, df = 1, log = TRUE,
     type = "shifted")
}
\arguments{
 \item{x}{Vector or matrix of quantiles. If \code{x} is a matrix, each
   row is taken to be a quantile.}
 \item{n}{Number of observations.}
 \item{delta}{ the vector of noncentrality parameters of length n, for
   \code{type = "shifted"} delta specifies the mode.}
 \item{sigma}{Correlation matrix, default is \code{diag(ncol(x))}.}
 \item{df}{ degree of freedom as integer.}
 \item{log}{Logical; if \code{TRUE}, densities d are given as log(d).}
 \item{type}{ type of the noncentral multivariate t distribution
              to be computed. \code{type = "Kshirsagar"} corresponds
              to formula (1.4) in Genz and Bretz (2009) (see also
	      Chapter 5.1 in Kotz and Nadarajah (2004)). This is the
	      noncentral t-distribution  needed for calculating
	      the power of multiple contrast tests under a normality
	      assumption. \code{type = "shifted"} corresponds to the
	      formula right before formula (1.4) in Genz and Bretz (2009)
	      (see also formula (1.1) in Kotz and Nadarajah (2004)). It
	      is a location shifted version of the central t-distribution.
	      This noncentral multivariate t distribution appears for
	      example as the Bayesian posterior distribution
	      for the regression coefficients in a linear regression.
 	      In the central case both types coincide.
              Note that the defaults differ from the default 
              in \code{\link{pmvt}} (for reasons of backward 
              compatibility).
	    }
 \item{\dots}{ additional arguments to \code{\link{rmvnorm}},
               for example \code{method}.}
	  }
\details{
  For type = "shifted" the following density is implemented
  \deqn{c(1+(x-\delta)'S^{-1}(x-\delta)/\nu)^{-(\nu+m)/2},}
  where
  \deqn{c = \Gamma((\nu+m)/2)/((\pi \nu)^{m/2}\Gamma(\nu/2)|S|^{1/2}),}
  here \eqn{S} is a positive definite symmetric matrix (which might be
  the correlation or the covariance matrix), \eqn{delta} is the
  non-centrality vector and  \eqn{\nu} are the degrees of freedom.
}
  
\seealso{\code{\link{pmvt}} and \code{\link{qmvt}}}
\examples{

  dmvt(x=c(0,0), sigma = diag(2))
  x <- rmvt(n=100, sigma = diag(2), df = 3)
  plot(x)

}
\keyword{distribution}
\keyword{multivariate}
