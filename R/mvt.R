# $Id: mvt.R,v 1.20 2002/04/09 14:16:48 hothorn Exp $ 

checkmvArgs <- function(lower, upper, mean, corr, sigma) 
{
    UNI <- FALSE
    if (is.null(lower) || any(is.na(lower)))
        stop("lower not specified or contains NA")
    if (is.null(upper) || any(is.na(upper)))
        stop("upper not specified or contains NA")
    rec <- cbind(lower, upper, mean)
    lower <- rec[,"lower"]
    upper <- rec[,"upper"]
    mean <- rec[,"mean"]
    if (any(is.na(mean)))
        stop("mean contains NA")
    if (is.null(corr) && is.null(sigma)) {
        corr <- diag(length(lower))
        # warning("both corr and sigma not specified: using sigma=diag(length(lower))")
    }
    if (!is.null(corr) && !is.null(sigma)) {
        sigma <- NULL
        warning("both corr and sigma specified: ignoring sigma")
    }
    if (!is.null(corr)) {
         if (!is.matrix(corr)) {
             if (length(corr) == 1)
                UNI <- TRUE
             if (length(corr) != length(lower))
               stop("diag(corr) and lower are of different length")
         } else {
             if (length(corr) == 1) {
                 UNI <- TRUE
                 corr <- corr[1,1]
             } else {
                 if (length(diag(corr)) != length(lower))        
                     stop("diag(corr) and lower are of different length")
             }
         }
    }
    if (!is.null(sigma)) {
         if (!is.matrix(sigma)) {
            if (length(sigma) == 1)
                UNI <- TRUE
            if (length(sigma) != length(lower))        
               stop("diag(sigma) and lower are of different length")
         } else {
            if (length(sigma) == 1) {
                UNI <- TRUE       
                sigma <- sigma[1,1]
            }
            if (length(diag(sigma)) != length(lower))                     
               stop("diag(sigma) and lower are of different length")
         }
    }
    list(lower=lower, upper=upper, mean=mean, corr=corr, sigma=sigma, uni=UNI)
}


pmvnorm <- function(lower=-Inf, upper=Inf, mean=rep(0, length(lower)), corr=NULL, sigma=NULL,
                    maxpts = 25000, abseps = 0.001, releps = 0)
{
    carg <- checkmvArgs(lower=lower, upper=upper, mean=mean, corr=corr,
                      sigma=sigma)
    if (!is.null(carg$corr)) {
      corr <- carg$corr
      if (carg$uni) {
          stop("sigma not specified: cannot compute pnorm")
      } else {
          lower <- carg$lower - carg$mean
          upper <- carg$upper - carg$mean
          mean <- rep(0, length(lower))
          RET <- mvt(lower=lower, upper=upper, df=0, corr=corr, delta=mean,
                     maxpts=maxpts, abseps=abseps,releps=releps)
      }
    } else {
      if (carg$uni) {
        RET <- list(value = pnorm(carg$upper, mean=carg$mean, sd=sqrt(carg$sigma)) -
                            pnorm(carg$lower, mean=carg$mean, sd=sqrt(carg$sigma)),
                    error = 0, msg="univariate: using pnorm")
      } else {
          lower <- (carg$lower - carg$mean)/sqrt(diag(carg$sigma))
          upper <- (carg$upper - carg$mean)/sqrt(diag(carg$sigma))
          mean <- rep(0, length(lower))
          corr <- sig2corr(carg$sigma)
          RET <- mvt(lower=lower, upper=upper, df=0, corr=corr, delta=mean,
                     maxpts=maxpts, abseps=abseps,releps=releps)
      }
    }
    attr(RET$value, "error") <- RET$error
    attr(RET$value, "msg") <- RET$msg
    return(RET$value)
}

pmvt <- function(lower=-Inf, upper=Inf, delta=rep(0, length(lower)),
                 df=1, corr=NULL, sigma=NULL, maxpts = 25000, abseps = 0.001,
                 releps = 0)
{
    carg <- checkmvArgs(lower=lower, upper=upper, mean=delta, corr=corr,
                       sigma=sigma)
    if (is.null(df))
        stop("df not specified")
    if (df < 1)
        stop("cannot compute multivariate t distribution with df < 1")
    if (carg$uni) {
        RET <- list(value = pt(carg$upper, df=df, ncp=carg$mean) -
                            pt(carg$lower, df=df, ncp=carg$mean),
                    error = 0, msg="univariate: using pt")
    } else {
        if (!is.null(carg$corr)) {
            RET <- mvt(lower=carg$lower, upper=carg$upper, df=df, corr=carg$corr,
                       delta=carg$mean,  maxpts=maxpts,
                       abseps=abseps,releps=releps)
        } else {
            lower <- carg$lower/sqrt(diag(carg$sigma))
            upper <- carg$upper/sqrt(diag(carg$sigma))
            corr <- sig2corr(carg$sigma)
            RET <- mvt(lower=lower, upper=upper, df=df, corr=corr,
                       delta=carg$mean, maxpts=maxpts,
                       abseps=abseps,releps=releps)
        }
    }
    attr(RET$value, "error") <- RET$error
    attr(RET$value, "msg") <- RET$msg
    return(RET$value)
}


mvt <- function(lower, upper, df, corr, delta, maxpts = 25000,
                abseps = 0.001, releps = 0)
{
    n <- ncol(corr)
    if (is.null(n) || n < 2) stop("dimension less then n = 2")

    if (length(lower) != n) stop("wrong dimensions")
    if (length(upper) != n) stop("wrong dimensions")

    if (n > 1000) stop("only dimensions 1 <= n <= 100 allowed") 

    infin <- rep(2, n)
    infin[upper == Inf] <- 1
    infin[lower == -Inf] <- 0
    infin[lower == -Inf & upper == Inf] <- -1
    
    if (n > 1) {
        corrF <- matrix(as.vector(corr), ncol=n, byrow=TRUE)
        corrF <- corrF[upper.tri(corrF)]
    } else corrF <- corr 

    lower[lower == -Inf] <- 0
    upper[upper == Inf] <- 0

    error <- 0; value <- 0; inform <- 0

    ret <- .Fortran("mvtdst", as.integer(n), as.integer(df),
                        as.double(lower), as.double(upper), as.integer(infin),
                        as.double(corrF), as.double(delta), as.integer(maxpts),
                        as.double(abseps), as.double(releps),  
                        error = as.double(error), value = as.double(value),
                        inform = as.integer(inform))
    
    error <- ret$error; value <- ret$value; inform <- ret$inform

    msg <- NULL
    if (inform == 0) msg <- "Normal Completion"
    if (inform == 1) msg <- "Completion with error > abseps"
    if (inform == 3) msg <- "Covariance matrix not positive semidefinite"
    if (is.null(msg)) msg <- inform
    
    RET <- list(value = value, error = error, msg = msg)
    return(RET)
}

