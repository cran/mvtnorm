
library("mvtnorm")

set.seed(29)

chk <- function(...) stopifnot(isTRUE(all.equal(..., check.attributes = FALSE)))

### N samples with N different covariance matrices

thischeck <- expression({
  N <- 10
  J <- 5
  lt <- ltMatrices(matrix(runif(N * J * (J + c(-1, 1)[dg + 1L]) / 2) + 1, ncol = N), 
                   diag = dg)
  lt <- ltMatrices(lt, diag = dg, byrow = br)
  Z <- matrix(rnorm(N * J), ncol = N)
  Y <- solve(lt, Z)
  ll1 <- sum(dnorm(Mult(lt, Y), log = TRUE)) + sum(log(diagonals(lt)))

  S <- as.array(Tcrossprod(solve(lt)))
  ll2 <- sum(l2 <- sapply(1:N, function(i) mvtnorm:::dmvnorm(x = Y[,i], sigma = S[,,i], log = TRUE)))
  chk(ll1, ll2)

  l3 <- ldmvnorm(obs = Y, invchol = lt, logLik = FALSE)
  l4 <- ldmvnorm(obs = Y, chol = solve(lt), logLik = FALSE)

  chk(l2, l3)
  chk(l2, l4)

  ll1 <- sum(dnorm(Mult(lt[1,], Y), log = TRUE)) + N * sum(log(diagonals(lt[1,])))

  S <- as.array(Tcrossprod(solve(lt)))
  ll2 <- sum(l2 <- sapply(1:N, function(i) mvtnorm:::dmvnorm(x = Y[,i], sigma = S[,,1], log = TRUE)))
  chk(ll1, ll2)

  l3 <- ldmvnorm(obs = Y, invchol = lt[1,], logLik = FALSE)
  l4 <- ldmvnorm(obs = Y, chol = solve(lt[1,]), logLik = FALSE)

  chk(l2, l3)
  chk(l2, l4)

  ### check scores
  if (require("numDeriv", quietly = TRUE)) {

    f <- function(L) {
      L <- ltMatrices(L, diag = dg, byrow = br)
      ldmvnorm(obs = Y, invchol = L)
    }

    s0 <- grad(f, unclass(lt))
    s1 <- sldmvnorm(obs = Y, invchol = lt)

    chk(Lower_tri(ltMatrices(matrix(s0, ncol = N), diag = dg, byrow = br), diag = dg), 
        Lower_tri(s1$invchol, diag = dg))

    f <- function(L) {
      L <- ltMatrices(L, diag = dg, byrow = br)
      ldmvnorm(obs = Y, chol = L)
    }

    s0 <- grad(f, unclass(lt))
    s1 <- sldmvnorm(obs = Y, chol = lt)

    chk(Lower_tri(ltMatrices(matrix(s0, ncol = N), diag = dg, byrow = br), diag = dg), 
        Lower_tri(s1$chol, diag = dg))

    f <- function(x)
      ldmvnorm(obs = x, invchol = lt)

    s0 <- grad(f, Y)
    s1 <- sldmvnorm(obs = Y, invchol = lt)

    chk(matrix(s0, ncol = N), s1$obs)
  }
})

dg <- TRUE
br <- FALSE
eval(thischeck)

dg <- FALSE
br <- FALSE
eval(thischeck)

dg <- FALSE
br <- TRUE
eval(thischeck)

dg <- FALSE
br <- FALSE
eval(thischeck)

