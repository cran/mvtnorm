
pmvt <- function(lower, upper, df, corr, delta, maxpts = 25000,
                 abseps = 0.001, releps = 0)
{
    if (df < 1) 
        stop("cannot compute multivariate t distribution with df < 1")
    if (is.null(ncol(corr)) || length(lower) == 1) 
        return(list(value = pt(upper, df=df, ncp=delta) -
                            pt(lower, df=df, ncp=delta),
                    error = 0, msg="univariate: using pt"))
    return(mvt(lower, upper, df, corr, delta, maxpts, abseps,releps))
}

pmvnorm <- function(lower, upper, mean, corr, maxpts = 25000,
                    abseps = 0.001, releps = 0)
{
    if (length(mean) == 2) {
      delta <- c(0,0)	# bug in the Fortran Sources FIXME!
      lower <- lower - mean
      upper <- upper - mean
    } else
      delta <- mean 
    if (length(mean) != length(lower)) stop("wrong dimensions")
    if (is.null(ncol(corr)) || length(mean) == 1)
        return(list(value = pnorm(upper, mean=mean, sd=corr) -
                            pnorm(lower, mean=mean, sd=corr),
                    error = 0, msg="univariate: using pnorm"))
    return(mvt(lower, upper, df=0, corr, delta, maxpts, abseps,releps))
}

mvt <- function(lower, upper, df, corr, delta, maxpts = 25000,
                abseps = 0.001, releps = 0)
{
    n <- ncol(corr)
    if (is.null(n) || n < 2) stop("dimension less then n = 2")

    if (length(lower) != n) stop("wrong dimensions")
    if (length(upper) != n) stop("wrong dimensions")

    if (n > 100) stop("only dimensions 1 <= n <= 100 allowed") 

    infin <- rep(2, n)
    infin[upper == Inf] <- 1
    infin[lower == -Inf] <- 0
    infin[lower == -Inf & upper == Inf] <- -1
    
    if (n > 1) {
        corrF <- matrix(as.vector(corr), ncol=n, byrow=T)
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
    
    out <- list(value = value, error = error, msg = msg)
    return(out)
}
