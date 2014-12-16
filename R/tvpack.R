
TVPACK <- function(abseps = 1e-6) {
    ret <- list(eps = abseps)
    class(ret) <- "TVPACK"
    ret
}

probval.TVPACK <- function (x, n, df, lower, upper, infin, corr, corrF, delta) {

    if (n > 3) 
      stop("TVPACK algorithms cannot compute probabilities for n > 3")
    if (df > 0 & any(delta != 0))
      stop("TVPACK only possible for the central t-distribution.")
  
    upp <- upper - delta
    low <- lower - delta
  
    if ((any(infin < 0) | any(infin > 1)) | length(unique(infin)) > 1)
        stop("TVPACK either needs all(lower == -Inf) or all(upper == Inf).")

    if (all(infin == 1))
        upp <- -low

    upp <- as.double(upp)
    eps <- as.double(x$eps)
    nu <- ifelse(df == 0, as.integer(0), as.integer(df))
  
    if (n == 2) {
        cr <- as.double(corr[2,1])
        res <- .C("C_bvtlr", nu, upp[1], upp[2], cr, val = double(1))
    }

    if (n == 3) {
        cr <- c(corr[2,1], corr[3,1], corr[3,2])
        cr <- as.double(cr)
        res <- .C("C_tvtlr", nu, upp, cr, eps, val = double(1))
    }
    error <- ifelse(n == 3, x$eps, NA)
    ret <- list(value = res$val, inform = 0, error = error)
    ret
}
