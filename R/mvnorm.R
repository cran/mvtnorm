# $Id: mvnorm.R,v 1.4 2004/01/22 14:30:50 hothorn Exp $

rmvnorm <- function(n, mean=rep(0, nrow(sigma)),
                      sigma=diag(length(mean))){

    if(nrow(sigma) != ncol(sigma)){
        stop("sigma must be a square matrix")
    }

    if(length(mean) != nrow(sigma)){
        stop("mean and sigma have non-conforming size")
    }

    ev <- eigen(sigma, sym = TRUE)$values
    if (!all(ev >= -sqrt(.Machine$double.eps) * abs(ev[1])))   
        warning("sigma is numerically not positive definite")

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
    logdet <- sum(log(eigen(sigma, symmetric=TRUE,
                                   only.values=TRUE)$values))
    logretval <- -(ncol(x)*log(2*pi) + logdet + distval)/2
    if(log) return(logretval)
    exp(logretval)
}
  
