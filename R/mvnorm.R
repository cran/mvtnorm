# $Id: mvnorm.R,v 1.2 2003/02/13 10:24:02 hothorn Exp $

rmvnorm <- function(n, mean=rep(0, nrow(sigma)),
                      sigma=diag(length(mean))){

  if(nrow(sigma) != ncol(sigma)){
    stop("sigma meanst be a square matrix")
  }

  if(length(mean) != nrow(sigma)){
    stop("mean and sigma have non-conforming size")
  }
  
  sigsvd <- svd(sigma)
  retval <- t(sigsvd$v %*% (t(sigsvd$u) * sqrt(sigsvd$d)))
  retval <- matrix(rnorm(n * ncol(sigma)), nrow = n) %*% retval
  retval <- sweep(retval, 2, mean, "+")
  retval
}


dmvnorm <- function (x, mean, sigma, log=FALSE)
{
    if (is.vector(x)) {
        x <- matrix(x, ncol = length(x))
    }
    if (missing(mean)) {
        mean <- rep(0, length = ncol(x))
    }
    if (missing(sigma)) {
        sigma <- diag(ncol(x))
    }
    if (ncol(x) != ncol(sigma)) {
        stop("x and sigma have non-conforming size")
    }
    if (nrow(sigma) != ncol(sigma)) {
        stop("sigma meanst be a square matrix")
    }
    if (length(mean) != nrow(sigma)) {
        stop("mean and sigma have non-conforming size")
    }
    distval <- mahalanobis(x, center = mean, cov = sigma)
    logdet <- sum(log(La.eigen(sigma, symmetric=TRUE,
                                      only.values=TRUE)$values))
    logretval <- -(ncol(x)*log(2*pi) + logdet + distval)/2
    if(log) return(logretval)
    exp(logretval)
}
  
