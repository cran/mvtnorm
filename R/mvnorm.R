# $Id: mvnorm.R 264 2014-04-04 10:19:18Z thothorn $

rmvnorm<-function (n, mean = rep(0, nrow(sigma)), sigma = diag(length(mean)),
                   method=c("eigen", "svd", "chol"), pre0.9_9994 = FALSE)
{    
    if (!isSymmetric(sigma, tol = sqrt(.Machine$double.eps), 
                     check.attributes = FALSE)) {
        stop("sigma must be a symmetric matrix")
    }
    if (length(mean) != nrow(sigma)) {
        stop("mean and sigma have non-conforming size")
    }
    sigma1 <- sigma
    dimnames(sigma1) <- NULL
    if(!isTRUE(all.equal(sigma1, t(sigma1)))){
        warning("sigma is numerically not symmetric")
    }

    method <- match.arg(method)
    
    if(method == "eigen"){
        ev <- eigen(sigma, symmetric = TRUE)
        if (!all(ev$values >= -sqrt(.Machine$double.eps) * abs(ev$values[1]))){
            warning("sigma is numerically not positive definite")
        }    
        retval <- ev$vectors %*%  diag(sqrt(ev$values), 
                      length(ev$values)) %*% t(ev$vectors)
    }
    else if(method == "svd"){
        sigsvd <- svd(sigma)
        if (!all(sigsvd$d >= -sqrt(.Machine$double.eps) * abs(sigsvd$d[1]))){
            warning("sigma is numerically not positive definite")
        }    
        retval <- t(sigsvd$v %*% (t(sigsvd$u) * sqrt(sigsvd$d)))
    }    
    else if(method == "chol"){
        retval <- chol(sigma, pivot = TRUE)
        o <- order(attr(retval, "pivot"))
        retval <- retval[,o]
    }
    
    retval <- matrix(rnorm(n * ncol(sigma)), nrow = n, byrow = !pre0.9_9994) %*%  retval
    retval <- sweep(retval, 2, mean, "+")
    colnames(retval) <- names(mean)
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
    if (NCOL(x) != NCOL(sigma)) {
        stop("x and sigma have non-conforming size")
    }
    if (!isSymmetric(sigma, tol = sqrt(.Machine$double.eps), 
                     check.attributes = FALSE)) {
        stop("sigma must be a symmetric matrix")
    }
    if (length(mean) != NROW(sigma)) {
        stop("mean and sigma have non-conforming size")
    }
    if( !is.null(dim(mean)) ) dim(mean) <- NULL

    ### <faster code contributed by Matteo Fasiolo mf364 at bath.ac.uk
    dec <- chol(sigma)
    tmp <- forwardsolve(t(dec), t(x) - mean)
    rss <- colSums(tmp ^ 2)
    logretval <- - sum(log(diag(dec))) - 0.5 * length(mean) * log(2 * pi) - 0.5 * rss
    ### />

    names(logretval) <- rownames(x)

    if(log) return(logretval)
    exp(logretval)
}
  