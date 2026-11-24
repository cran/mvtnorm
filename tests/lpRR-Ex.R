
library("mvtnorm")

set.seed(29081975)

chk <- function(...)
    stopifnot(all.equal(..., tol = 1e-5, check.attributes = FALSE))

### lpRR (reduced rank, Sigma = B B' + D) and lpmvnorm (Genz, full rank on
### the same Sigma) integrate the same probability two different ways -- they
### should agree, and so should their scores.  tests/regtest-lpRR.R covers
### the numerical edge cases (underflow, missing limits, tol); this file is
### the minimal cross-implementation demonstration, also far in the tails.

J <- 6
K <- 3
B <- matrix(rnorm(J * K), nrow = J) / sqrt(K)
D <- runif(J, .3, 1)
mn <- rnorm(J) / 4
a <- -(1 + runif(J))
b <- 1 + runif(J)

Sigma <- tcrossprod(B) + diag(D)
C <- t(chol(Sigma))
Cl <- as.chol(ltMatrices(C[lower.tri(C, diag = TRUE)], diag = TRUE))

M <- 1e5
Z <- matrix(rnorm(K * M), nrow = K)
w <- matrix(runif((J - 1) * M), nrow = J - 1, byrow = TRUE)

aM <- matrix(a, ncol = 1); bM <- matrix(b, ncol = 1); mnM <- matrix(mn, ncol = 1)

### log-probabilities agree, up to Monte Carlo error at this M

lRR <- lpRR(lower = a, upper = b, mean = mn, B = B, D = D, Z = Z)
lGz <- lpmvnorm(lower = aM, upper = bM, mean = mnM, chol = Cl, w = w, M = M)
chk(lRR, lGz, tolerance = 1e-2)

### scores agree on the parameters the two parameterisations share: lower,
### upper, mean.  B/D vs chol are different parameterisations of the same
### Sigma and are not directly comparable.

sRR <- slpRR(lower = a, upper = b, mean = mn, B = B, D = D, Z = Z)
sGz <- slpmvnorm(lower = aM, upper = bM, mean = mnM, chol = Cl, w = w, M = M)
chk(sRR$lower, c(sGz$lower), tolerance = 1e-2)
chk(sRR$upper, c(sGz$upper), tolerance = 1e-2)
chk(sRR$mean,  c(sGz$mean),  tolerance = 1e-2 * 2) # mean is small in magnitude,
                                                    # so relative tolerance needs slack

### and both scores actually belong to their own log-likelihood: analytic
### score vs numDeriv::grad(), with the SAME integration points (Z or w)
### inside the numeric and the analytic call

### speed up, as precision of log-prob is not relevant for evaluating
### the gradient function
M <- 1e3
Z <- matrix(rnorm(K * M), nrow = K)
sRR <- slpRR(lower = a, upper = b, mean = mn, B = B, D = D, Z = Z)

if (require("numDeriv", quietly = TRUE)) {

    llRR <- function(p, what) {
        A <- a; U <- b; Mu <- mn; BB <- B; DD <- D
        switch(what, lower = A <- p, upper = U <- p, mean = Mu <- p,
               B = BB <- matrix(p, nrow = J), D = DD <- p)
        lpRR(lower = A, upper = U, mean = Mu, B = BB, D = DD, Z = Z)
    }
    chk(grad(llRR, a, what = "lower"), c(sRR$lower))
    chk(grad(llRR, b, what = "upper"), c(sRR$upper))
    chk(grad(llRR, mn, what = "mean"), c(sRR$mean))
    chk(grad(llRR, c(B), what = "B"), c(sRR$B))
    chk(grad(llRR, D, what = "D"), c(sRR$D))

}

### near-degenerate D: lpRR/slpRR stay finite -- the whole point of this
### rewrite (see lmvnorm_src.w, chapter "Reduced Rank Covariance Matrices").

Dd <- c(rep(1e-3, 3), D[-(1:3)])
### 1e5 points as in w, for the comparison with lpmvnorm below
Zd <- matrix(rnorm(K * ncol(w)), nrow = K)
lRRd <- lpRR(lower = a, upper = b, mean = mn, B = B, D = Dd, Z = Zd)
stopifnot(is.finite(lRRd))
sRRd <- slpRR(lower = a, upper = b, mean = mn, B = B, D = Dd, Z = Zd)
stopifnot(all(sapply(sRRd, function(x) all(is.finite(x)))))

### lpmvnorm/slpmvnorm agree with lpRR/slpRR in this regime, too

Sigmad <- tcrossprod(B) + diag(Dd)
Cd <- t(chol(Sigmad))
Cld <- as.chol(ltMatrices(Cd[lower.tri(Cd, diag = TRUE)], diag = TRUE))
lGzd <- lpmvnorm(lower = aM, upper = bM, mean = mnM, chol = Cld, w = w)
chk(lRRd, lGzd, tolerance = 1e-1)
sGzd <- slpmvnorm(lower = aM, upper = bM, mean = mnM, chol = Cld, w = w)
chk(sRRd$lower, c(sGzd$lower), tolerance = 1e-1)
chk(sRRd$upper, c(sGzd$upper), tolerance = 1e-1)
chk(sRRd$mean,  c(sGzd$mean),  tolerance = 1e-1)

### far in the tail: lpmvnorm vs lpRR with K = 1 on a fine grid (exact
### quadrature); boxes whose signs follow the loadings are avoided (large
### lattice error, grid weights underflow)

h <- .005
grid <- seq(-60, 60, by = h)
Zg <- matrix(grid, nrow = 1)
wg <- dnorm(grid) * h

J1 <- 5
B1 <- matrix(c(.8, -.6, .5, -.7, .4), ncol = 1)
D1 <- c(.5, .8, .3, .6, .9)
C1 <- t(chol(tcrossprod(B1) + diag(D1)))
C1 <- ltMatrices(C1[lower.tri(C1, diag = TRUE)], diag = TRUE)
mn1 <- c(.2, -.1, 0, .3, -.2)
M1 <- 20000
w1 <- (outer(sqrt(c(2, 3, 5, 7)), 1:M1) + runif(J1 - 1)) %% 1

for (shift in c(-30, -10, 10, 30)) {
    a1 <- shift - c(.2, .5, .8, .1, .6)
    b1 <- a1 + c(1.5, 1, 2, 1.2, 1.8)
    r <- lpRR(a1, b1, mn1, B = B1, D = D1, Z = Zg, weights = wg)
    n <- lpmvnorm(matrix(a1), matrix(b1), mean = matrix(mn1), chol = C1, w = w1)
    sr <- slpRR(a1, b1, mn1, B = B1, D = D1, Z = Zg, weights = wg)
    sn <- slpmvnorm(matrix(a1), matrix(b1), mean = matrix(mn1), chol = C1, w = w1)
    stopifnot(is.finite(r), abs(n - r) < .01)
    nms <- c("lower", "upper", "mean")
    rel <- mean(abs(unlist(sr[nms]) - unlist(lapply(sn[nms], c)))) /
           mean(abs(unlist(sr[nms])))
    stopifnot(rel < 1e-4)
}
