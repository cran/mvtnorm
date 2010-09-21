# $Id: mvt.R 219 2010-11-16 17:18:01Z thothorn $ 

chkcorr <- function(x) {

    if (!is.matrix(x)) return(TRUE)
    rownames(x) <- colnames(x) <- NULL
    storage.mode(x) <- "numeric"
    ONE <- 1 + .Machine$double.eps

    ret <- (min(x) < -ONE || max(x) > ONE) ||
           !isTRUE(all.equal(diag(x), rep(1, nrow(x))))
    !ret
}

checkmvArgs <- function(lower, upper, mean, corr, sigma) 
{
    UNI <- FALSE
    if (!is.numeric(lower) || !is.vector(lower))
        stop(sQuote("lower"), " is not a numeric vector")
    if (!is.numeric(upper) || !is.vector(upper))
        stop(sQuote("upper"), " is not a numeric vector")
    if (!is.numeric(mean) || !is.vector(mean))
        stop(sQuote("mean"), " is not a numeric vector")
    if (is.null(lower) || any(is.na(lower)))
        stop(sQuote("lower"), " not specified or contains NA")
    if (is.null(upper) || any(is.na(upper)))
        stop(sQuote("upper"), " not specified or contains NA")
    rec <- cbind(lower, upper, mean)
    lower <- rec[,"lower"]
    upper <- rec[,"upper"]
    if (!all(lower <= upper))
        stop("at least one element of ", sQuote("lower"), " is larger than ", 
             sQuote("upper"))
    mean <- rec[,"mean"]
    if (any(is.na(mean)))
        stop("mean contains NA")
    if (is.null(corr) && is.null(sigma)) {
        corr <- diag(length(lower))
        # warning("both ", sQuote("corr"), " and ", sQuote("sigma"),
        # " not specified: using sigma=diag(length(lower))")
    }
    if (!is.null(corr) && !is.null(sigma)) {
        sigma <- NULL
        warning("both ", sQuote("corr"), " and ", sQuote("sigma"), 
                " specified: ignoring ", sQuote("sigma"))
    }
    if (!is.null(corr)) {
         if (!is.numeric(corr))
             stop(sQuote("corr"), " is not numeric")
         if (!is.matrix(corr)) {
             if (length(corr) == 1)
                UNI <- TRUE
             if (length(corr) != length(lower))
               stop(sQuote("diag(corr)"), " and ", sQuote("lower"), 
                    " are of different length")
         } else {
             if (length(corr) == 1) {
                 UNI <- TRUE
                 corr <- corr[1,1]
                 if (length(lower) != 1)
                   stop(sQuote("corr"), " and ", sQuote("lower"), 
                        " are of different length")
             } else {
                 if (length(diag(corr)) != length(lower))        
                     stop(sQuote("diag(corr)"), " and ", sQuote("lower"), 
                          " are of different length")
                 if (!chkcorr(corr))
                     stop(sQuote("corr"), " is not a correlation matrix")
             }
         }
    }
    if (!is.null(sigma)) {
         if (!is.numeric(sigma))
             stop(sQuote("sigma"), " is not numeric")
         if (!is.matrix(sigma)) {
            if (length(sigma) == 1)
                UNI <- TRUE
            if (length(sigma) != length(lower))        
               stop(sQuote("diag(sigma)"), " and ", sQuote("lower"), 
                    " are of different length")
         } else {
            if (length(sigma) == 1) {
                UNI <- TRUE       
                sigma <- sigma[1,1]
                if (length(lower) != 1) 
                  stop(sQuote("sigma"), " and ", sQuote("lower"), 
                       " are of different length")
            } else {
              if (length(diag(sigma)) != length(lower))                     
                 stop(sQuote("diag(sigma)"), " and ", sQuote("lower"), 
                      " are of different length")
              if (!isTRUE(all.equal(sigma, t(sigma))) || any(diag(sigma) < 0))
                 stop(sQuote("sigma"), " is not a covariance matrix")
            }
         }
    }
    list(lower=lower, upper=upper, mean=mean, corr=corr, sigma=sigma, uni=UNI)
}


pmvnorm <- function(lower=-Inf, upper=Inf, mean=rep(0, length(lower)), corr=NULL, sigma=NULL,
                    algorithm = GenzBretz(), ...)
{
    carg <- checkmvArgs(lower=lower, upper=upper, mean=mean, corr=corr,
                        sigma=sigma)
    if (!is.null(carg$corr)) {
      corr <- carg$corr
      if (carg$uni) {
          stop(sQuote("sigma"), " not specified: cannot compute pnorm")
      } else {
          lower <- carg$lower - carg$mean
          upper <- carg$upper - carg$mean
          mean <- rep(0, length(lower))
          RET <- mvt(lower=lower, upper=upper, df=0, corr=corr, delta=mean,
                     algorithm = algorithm, ...)
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
          corr <- cov2cor(carg$sigma)
          RET <- mvt(lower=lower, upper=upper, df=0, corr=corr, delta=mean,
                     algorithm = algorithm, ...)
      }
    }
    attr(RET$value, "error") <- RET$error
    attr(RET$value, "msg") <- RET$msg
    return(RET$value)
}

pmvt <- function(lower=-Inf, upper=Inf, delta=rep(0, length(lower)),
                 df=1, corr=NULL, sigma=NULL, 
                 algorithm = GenzBretz(), 
                 type = c("Kshirsagar", "shifted"), ...)
{
    type <- match.arg(type)
    carg <- checkmvArgs(lower=lower, upper=upper, mean=delta, corr=corr,
                        sigma=sigma)
    if (type == "shifted") { # can be handled by integrating over central t
      if(!is.null(carg$corr)){ # using transformed integration bounds
        d <- 1
      } else {
        if(!is.null(carg$sigma)){
          d <- sqrt(diag(carg$sigma))
          carg$corr <- cov2cor(carg$sigma)
        }
      }
      carg$lower <- (carg$lower - carg$mean)/d
      carg$upper <- (carg$upper - carg$mean)/d
      carg$mean <- rep(0, length(carg$mean))
    }

    if (is.null(df))
        stop(sQuote("df"), " not specified")
    if (any(df < 0))
        stop("cannot compute multivariate t distribution with ",
             sQuote("df"), " < 0")
    if (carg$uni) {
        if (df > 0)
            RET <- list(value = pt(carg$upper, df=df, ncp=carg$mean) -
                                pt(carg$lower, df=df, ncp=carg$mean),
                       error = 0, msg="univariate: using pt")
        else
            RET <- list(value = pnorm(carg$upper, mean = carg$mean) -
                                pnorm(carg$lower, mean=carg$mean),
                       error = 0, msg="univariate: using pnorm")
    } else {
        if (!is.null(carg$corr)) {
            RET <- mvt(lower=carg$lower, upper=carg$upper, df=df, corr=carg$corr,
                       delta=carg$mean,  algorithm = algorithm, ...)
        } else { # need to transform integration bounds and delta
            d <- sqrt(diag(carg$sigma))
            lower <- carg$lower/d
            upper <- carg$upper/d
            corr <- cov2cor(carg$sigma)
            RET <- mvt(lower=lower, upper=upper, df=df, corr=corr,
                       delta=carg$mean/d, algorithm = algorithm, ...)
        }
    }
    attr(RET$value, "error") <- RET$error
    attr(RET$value, "msg") <- RET$msg
    return(RET$value)
}


mvt <- function(lower, upper, df, corr, delta, algorithm = GenzBretz(), ...)
{

    ### only for compatibility with older versions
    addargs <- list(...)
    if (length(addargs) > 0)
        algorithm <- GenzBretz(...)
    if (is.function(algorithm) || is.character(algorithm))
        algorithm <- do.call(algorithm, list())

    ### handle cases where the support is the empty set
    if (any(abs(lower - upper) < sqrt(.Machine$double.eps)) || 
        any(is.na(lower - upper))) {
        RET <- list(value = 0, error = 0, msg = "lower == upper")
        return(RET)
    }

    n <- ncol(corr)
    if (is.null(n) || n < 2) stop("dimension less then n = 2")

    if (length(lower) != n) stop("wrong dimensions")
    if (length(upper) != n) stop("wrong dimensions")

    if (n > 1000) stop("only dimensions 1 <= n <= 1000 allowed") 

    infin <- rep(2, n)   
    infin[upper == Inf] <- 1
    infin[lower == -Inf] <- 0
    infin[lower == -Inf & upper == Inf] <- -1

    ### this is a bug in `mvtdst' not yet fixed
    if (all(infin < 0)) 
        return(list(value = 1, error = 0, msg = "Normal Completion"))
    
    if (n > 1) {
        corrF <- matrix(as.vector(corr), ncol=n, byrow=TRUE)
        corrF <- corrF[upper.tri(corrF)]
    } else corrF <- corr 


    ret <- probval(algorithm, n, df, lower, upper, infin, corr, corrF, delta)
    error <- ret$error; value <- ret$value; inform <- ret$inform

    msg <- NULL
    if (inform == 0) msg <- "Normal Completion"
    if (inform == 1) msg <- "Completion with error > abseps"
    if (inform == 2) msg <- "N greater 1000 or N < 1"
    if (inform == 3) msg <- "Covariance matrix not positive semidefinite"
    if (is.null(msg)) msg <- inform
    
    RET <- list(value = value, error = error, msg = msg)
    return(RET)
}

rmvt <- function(n, sigma = diag(2), df = 1, 
     delta = rep(0, nrow(sigma)),
     type = c("shifted", "Kshirsagar")) {

    if (length(delta) != nrow(sigma))
      stop("delta and sigma have non-conforming size")

    if (df == 0)
        return(rmvnorm(n, mean = delta, sigma = sigma))
    type <- match.arg(type)

    if (type == "Kshirsagar")
        return(rmvnorm(n, mean = delta, sigma = sigma)/
               sqrt(rchisq(n, df)/df))

    if (type == "shifted"){
        sims <- rmvnorm(n, sigma = sigma)/sqrt(rchisq(n, df)/df)
        return(sweep(sims, 2, delta, "+"))
    }
    ### was: rmvnorm(n,sigma=sigma)/sqrt(rchisq(n,df)/df)
}

dmvt <- function(x, delta, sigma, df = 1, 
                 log = TRUE, type = "shifted")
{
    if (df == 0)
        return(dmvnorm(x, mean = delta, sigma = sigma, log = log))

    if (is.vector(x)) {
        x <- matrix(x, ncol = length(x))
    }
    if (missing(delta)) {
        delta <- rep(0, length = ncol(x))
    }
    if (missing(sigma)) {
        sigma <- diag(ncol(x))
    }
    if (NCOL(x) != NCOL(sigma)) {
        stop("x and sigma have non-conforming size")
    }
    if (!isSymmetric(sigma, tol = sqrt(.Machine$double.eps))) {
        stop("sigma must be a symmetric matrix")
    }
    if (length(delta) != NROW(sigma)) {
        stop("mean and sigma have non-conforming size")
    }

    m <- NCOL(sigma)

    distval <- mahalanobis(x, center = delta, cov = sigma)

    logdet <- sum(log(eigen(sigma, symmetric = TRUE,
                            only.values = TRUE)$values))

    logretval <- lgamma((m + df)/2) - 
                 (lgamma(df / 2) + 0.5 * (logdet + m * logb(pi * df))) -
                 0.5 * (df + m) * logb(1 + distval / df)
    if (log)
        return(logretval) 
    return(exp(logretval))
}

### find suitable interval to search for quantile
approx_interval <- function(p, tail, corr, df = 0) {

    qfun <- function(p) {
        if (df == 0) return(qnorm(p))
        return(qt(p, df = df))
    }

    if (tail == "both.tails")
        p <- p + (1 - p) / 2

             ### univariate quantile (corr == 1, perfect correlation)
    ret <- c(qfun(p),
             ### multivariate quantile (corr == 0, independence)
             qfun(p^(1 / NCOL(corr))))

    ### just to please uniroot
    ret <- ret * c(0.9, 1.1)

    if (tail == "upper.tail")
        ret <- rev(-ret)
    ret
}

qmvnorm <- function(p, interval = NULL, 
                    tail = c("lower.tail", "upper.tail", "both.tails"), 
                    mean = 0, corr = NULL, sigma = NULL, algorithm = 
                    GenzBretz(), ...)
{
    if (length(p) != 1 || (p <= 0 || p >= 1)) 
        stop(sQuote("p"), " is not a double between zero and one")

    dots <- dots2GenzBretz(...)
    if (!is.null(dots$algorithm) && !is.null(algorithm)) 
        algorithm <- dots$algorithm

    tail <- match.arg(tail)
    if (tail == "both.tails" && p < 0.5)
        stop("cannot compute two-sided quantile for p < 0.5")
    dim <- length(mean)
    if (is.matrix(corr)) dim <- nrow(corr)
    if (is.matrix(sigma)) dim <- nrow(sigma)
    lower <- rep(0, dim)
    upper <- rep(0, dim)
    args <- checkmvArgs(lower, upper, mean, corr, sigma)
    if (args$uni) {
        if (is.null(args$sigma))
          stop(sQuote("sigma"), " not specified: cannot compute qnorm")
        if (tail == "both.tails") p <- ifelse(p < 0.5, p / 2, 1 - (1 - p)/2)
        q <- qnorm(p, mean = args$mean, sd = args$sigma, 
                   lower.tail = (tail != "upper.tail"))
        qroot <- c(quantile = q, f.quantile = 0)
        return(qroot)
    }
    dim <- length(args$mean)

    pfct <- function(q) {
        switch(tail, "both.tails" = {
                  low <- rep(-abs(q), dim)
                  upp <- rep( abs(q), dim)
           }, "upper.tail" = {
                  low <- rep(      q, dim)
                  upp <- rep(    Inf, dim)
           }, "lower.tail" = {
                  low <- rep(   -Inf, dim)
                  upp <- rep(      q, dim)
           },)
           pmvnorm(lower = low, upper = upp, mean = args$mean,
                   corr = args$corr, sigma = args$sigma, 
                   algorithm = algorithm) - p
    }

    if (is.null(interval) || length(interval) != 2) {
        if (mean == 0) {
            cr <- args$corr
            if (!is.null(args$sigma)) cr <- cov2cor(args$sigma)
            interval <- approx_interval(p = p, tail = tail, 
                                        corr = cr, df = 0)
        } else {
            interval <- c(ifelse(tail == "both.tails", 0, -10), 10)
        }
    }

    if (is.null(dots$uniroot)) {
        qroot <- uniroot(pfct, interval = interval)
    } else {
        qroot <- do.call("uniroot", list(pfct, interval = interval, dots$uniroot))
    }
    names(qroot)[1:2] <- c("quantile", "f.quantile")
    qroot
}

qmvt <- function(p, interval = NULL, 
                 tail = c("lower.tail", "upper.tail", "both.tails"), 
                 df = 1, delta = 0, corr = NULL, sigma = NULL,
                 algorithm = GenzBretz(), 
                 type = c("Kshirsagar", "shifted"), ...) {

    if (length(p) != 1 || (p <= 0 || p >= 1)) 
        stop(sQuote("p"), " is not a double between zero and one")

    dots <- dots2GenzBretz(...)
    if (!is.null(dots$algorithm)  && !is.null(algorithm)) 
        algorithm <- dots$algorithm
    type <- match.arg(type)

    tail <- match.arg(tail)
    if (tail == "both.tails" && p < 0.5)
        stop("cannot compute two-sided quantile for p < 0.5")
    dim <- 1
    if (!is.null(corr)) dim <- NROW(corr)
    if (!is.null(sigma)) dim <- NROW(sigma)
    lower <- rep(0, dim)
    upper <- rep(0, dim)
    args <- checkmvArgs(lower, upper, delta, corr, sigma)
    if (args$uni) {
        if (tail == "both.tails") p <- ifelse(p < 0.5, p / 2, 1 - (1 - p)/2)
        q <- qt(p, df = df, ncp = args$mean, lower.tail = (tail != "upper.tail"))
        qroot <- c(quantile = q, f.quantile = 0)
        return(qroot)
    }

    dim <- length(args$mean)

    pfct <- function(q) {
        switch(tail, "both.tails" = {
                  low <- rep(-abs(q), dim)
                  upp <- rep( abs(q), dim)
           }, "upper.tail" = {
                  low <- rep(      q, dim)
                  upp <- rep(    Inf, dim)
           }, "lower.tail" = {
                  low <- rep(   -Inf, dim)
                  upp <- rep(      q, dim)
           },)
           pmvt(lower = low, upper = upp, df = df, delta = args$mean,
                corr = args$corr, sigma = args$sigma,
                algorithm = algorithm, type = type) - p
    }

    if (is.null(interval) || length(interval) != 2) {
        if (delta == 0) {
            cr <- args$corr
            if (!is.null(args$sigma)) cr <- cov2cor(args$sigma)
            interval <- approx_interval(p = p, tail = tail, 
                                        corr = cr, df = df)
        } else {
            interval <- c(ifelse(tail == "both.tails", 0, -10), 10)
        } 
    }

    if (is.null(dots$uniroot)) {
        qroot <- uniroot(pfct, interval = interval)
    } else {
        qroot <- do.call("uniroot", list(pfct, interval = interval, dots$uniroot))
    }
    names(qroot)[1:2] <- c("quantile", "f.quantile")
    qroot
}

GenzBretz <- function(maxpts = 25000, abseps = 0.001, releps = 0) {
    ret <- list(maxpts = maxpts, abseps = abseps, releps = releps)
    class(ret) <- "GenzBretz"
    ret
}

Miwa <- function(steps = 128) {
    ret <- list(steps = steps)
    class(ret) <- "Miwa"
    ret
}

probval <- function(x, ...)
    UseMethod("probval")

probval.GenzBretz <- function(x, n, df, lower, upper, infin, corr, corrF, delta) {

    lower[lower == -Inf] <- 0
    upper[upper == Inf] <- 0

    error <- 0; value <- 0; inform <- 0
    ret <- .Fortran("mvtdst", N = as.integer(n),
                              NU = as.integer(df),
                              LOWER = as.double(lower), 
                              UPPER = as.double(upper),
                              INFIN = as.integer(infin),
                              CORREL = as.double(corrF),
                              DELTA = as.double(delta),
                              MAXPTS = as.integer(x$maxpts),
                              ABSEPS = as.double(x$abseps),
                              RELEPS = as.double(x$releps),
                              error = as.double(error),
                              value = as.double(value),
                              inform = as.integer(inform), PACKAGE="mvtnorm")
    ret
}

probval.Miwa <- function(x, n, df, lower, upper, infin, corr, corrF, delta) {

    if (df != 0)
        stop("Miwa algorithm cannot compute t-probabilities")

    if (n > 20) 
        stop("Miwa algorithm cannot compute exact probabilities for n > 20")

    sc <- try(solve(corr))
    if (inherits(sc, "try-error")) 
        stop("Miwa algorithm cannot compute probabilities for singular problems")

    p <- .Call("C_miwa", steps = as.integer(x$steps),
                         corr = as.double(corr),
                         upper = as.double(upper),
                         lower = as.double(lower),
                         infin = as.integer(infin))
    ret <- list(value = p, inform = 0, error = NA)
    ret
}

dots2GenzBretz <- function(...) {
    addargs <- list(...)
    fm1 <- sapply(names(addargs), function(x) length(grep(x, names(formals(GenzBretz)))) == 1)
    fm2 <- sapply(names(addargs), function(x) length(grep(x, names(formals(uniroot)))) == 1)
    algorithm <- NULL
    uniroot <- NULL
    if (any(fm1))
        algorithm <- do.call("GenzBretz", addargs[fm1])
    if (any(fm2))
        uniroot <- addargs[fm2]
    list(algorithm = algorithm, uniroot = uniroot)
}
