sig2corr <- function(sigma) {
  if (!is.matrix(sigma)) stop("sigma is not a matrix")
  if (ncol(sigma) != nrow(sigma)) stop("sigma is not a square matrix")
  if (!all.equal(sigma, t(sigma))) stop("sigma is not symmetric")
  dummy <- matrix(1, ncol=ncol(sigma), nrow=nrow(sigma))
  dummy <- dummy/sqrt(diag(sigma))
  dummy <- t(t(dummy)/sqrt(diag(sigma)))
  corr <- dummy*sigma
  if (!all.equal(sum(diag(corr)), nrow(sigma))) stop("error in sig2corr")
  corr
}
