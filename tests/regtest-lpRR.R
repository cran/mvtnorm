
library("mvtnorm")

set.seed(29081975)

chk <- function(...)
    stopifnot(all.equal(..., tol = 1e-5, check.attributes = FALSE))

EVAL <- function(...) {}

if (require("numDeriv", quietly = TRUE))
  EVAL <- eval

J <- 6
K <- 3
B <- matrix(rnorm(J * K), nrow = J) / sqrt(K)
D <- runif(J, .3, 1)
S <- tcrossprod(B) + diag(D)
a <- -(1 + runif(J))
b <- 1 + runif(J)
M <- 1e5
Z <- matrix(rnorm(K * M), nrow = K)

ref <- function(lower, upper, Sigma)
    log(pmvnorm(lower = lower, upper = upper, sigma = Sigma,
                algorithm = GenzBretz(maxpts = 1e6, abseps = 0), seed = 1))

### independent oracle: lpRR vs pmvnorm on Sigma = B B' + D

stopifnot(abs((lp0 <- lpRR(lower = a, upper = b, B = B, D = D, Z = Z)) -
              ref(a, b, S)) < 1e-2)

### probabilities and log-probabilities agree

chk(exp(lp0),
    lpRR(lower = a, upper = b, B = B, D = D, Z = Z, log.p = FALSE))

### missing limits are the corresponding infinities
### (`upper' missing used to fail with "object 'pu' not found")

chk(lpRR(lower = a, B = B, D = D, Z = Z),
    lpRR(lower = a, upper = rep(Inf, J), B = B, D = D, Z = Z))
chk(lpRR(upper = b, B = B, D = D, Z = Z),
    lpRR(lower = rep(-Inf, J), upper = b, B = B, D = D, Z = Z))

### scores against their numerical approximations, with the SAME points
### inside the numeric and the analytic call

Zs <- matrix(rnorm(K * 2000), nrow = K)
mn <- rnorm(J) / 4

thischeck <- quote({

    sc <- slpRR(lower = a, upper = b, mean = mn, B = B, D = D, Z = Zs)

    ll <- function(p, what) {
        A <- a; U <- b; Mu <- mn; BB <- B; DD <- D
        switch(what, lower = A <- p, upper = U <- p, mean = Mu <- p,
               B = BB <- matrix(p, nrow = J), D = DD <- p)
        lpRR(lower = A, upper = U, mean = Mu, B = BB, D = DD, Z = Zs)
    }

    chk(grad(ll, a, what = "lower"), c(sc$lower))
    chk(grad(ll, b, what = "upper"), c(sc$upper))
    chk(grad(ll, mn, what = "mean"), c(sc$mean))
    chk(grad(ll, c(B), what = "B"), c(sc$B))
    chk(grad(ll, D, what = "D"), c(sc$D))
})

EVAL(thischeck)

### near-degenerate D: each factor pnorm(u) - pnorm(l) underflows to an exact
### zero on the linear scale, which used to annihilate the whole point and
### return -Inf.  The log-scale evaluation keeps these finite.

Dd <- c(rep(1e-3, 3), D[-(1:3)])
Sd <- tcrossprod(B) + diag(Dd)

### for the record: some cells of the old linear-scale integrant ARE zero
BZ <- (B / sqrt(Dd)) %*% Z
inner <- pmax(0, pnorm(b / sqrt(Dd) - BZ) - pnorm(a / sqrt(Dd) - BZ))
stopifnot(any(inner == 0))

pd <- lpRR(lower = a, upper = b, B = B, D = Dd, Z = Z)
stopifnot(is.finite(pd))
stopifnot(abs(pd - ref(a, b, Sd)) < 1e-1)

sd0 <- slpRR(lower = a, upper = b, B = B, D = Dd, Z = Z)
stopifnot(all(sapply(sd0, function(x) all(is.finite(x)))))

### no underflow floor: lpRR returns the true, correctly computed (very
### negative but finite) log-probability for a box far in the tail - not a
### floored approximation and not a cancellation artefact, since
### .lognormdiff already handles the per-coordinate tail on its own.  slpRR
### returns the correspondingly large, but finite, gradient.

far <- rep(60, J)
lfar <- lpRR(lower = far, upper = far + 1, B = B, D = D, Z = Z)
stopifnot(is.finite(lfar) && lfar < -100)
sfar <- slpRR(lower = far, upper = far + 1, B = B, D = D, Z = Z)
stopifnot(all(sapply(sfar, function(x) all(is.finite(x)))))
stopifnot(any(sapply(sfar, function(x) any(x != 0))))

### a degenerate box (both limits at the same infinity) is empty, so the
### probability is exactly zero; lpRR returns -Inf (this used to error out of
### .log1mexp with "NAs are not allowed in subscripted assignments").  No
### score is defined for a zero probability: slpRR errors rather than
### silently returning a zero score.

stopifnot(lpRR(lower = rep(Inf, J), upper = rep(Inf, J),
               B = B, D = D, Z = Z) == -Inf)
stopifnot(lpRR(lower = rep(-Inf, J), upper = rep(-Inf, J),
               B = B, D = D, Z = Z) == -Inf)
stopifnot(inherits(tryCatch(slpRR(lower = rep(Inf, J), upper = rep(Inf, J),
                                  B = B, D = D, Z = Z),
                            error = function(e) e), "error"))

