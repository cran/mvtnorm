### Tests for the log-scale implementation of lpmvnorm() and slpmvnorm()
### (1) exact values and scores far in the tails, independent case
### (2) scores vs numerical derivatives far in the tails, correlated case
### (3) degenerate situations
### (correlated case vs lpRR() far in the tails: lpRR-Ex.R)

library("mvtnorm")

set.seed(290875)

chk <- function(...) stopifnot(isTRUE(all.equal(...)))
allfinite <- function(x) all(sapply(x, function(x) all(is.finite(unclass(x)))))

### exact univariate oracle for N(mean, sd^2): log(P(l < Y < u)) computed
### stably in the tails (by reflection) and its scores. Returns a list of
### ll, lower, upper, mean (all vectors of the length of the limits)
uni <- function(l, u, m = 0, s = 1) {
    a <- (l - m) / s
    b <- (u - m) / s
    ### reflect the upper tail
    ra <- ifelse(a > 0, -b, a)
    rb <- ifelse(a > 0, -a, b)
    lb <- pnorm(rb, log.p = TRUE)
    la <- pnorm(ra, log.p = TRUE)
    ll <- ifelse(ra < 0 & rb > 0, log(pnorm(rb) - pnorm(ra)),
                 lb + log1p(-exp(la - lb)))
    dl <- -exp(dnorm(a, log = TRUE) - ll) / s
    du <- exp(dnorm(b, log = TRUE) - ll) / s
    list(ll = ll, lower = dl, upper = du, mean = -(dl + du))
}

### 1. exact values, independent case

### diagonal covariance: the exact answer is a sum of univariate terms
### (for any nodes w)
J <- 5
w <- matrix(runif((J - 1) * 10), nrow = J - 1)
for (shift in c(0, 10, 20, 35))
for (sg in list(rep(-1, J), rep(1, J), c(1, -1, 1, -1, 1))) {
    D <- runif(J, .3, 1)
    mn <- rnorm(J, sd = .3)
    lwr <- sg * shift - runif(J)
    upr <- lwr + runif(J, 1, 2)
    ref <- uni(lwr, upr, mn, sqrt(D))
    C <- ltMatrices(rep(0, J * (J + 1) / 2), diag = TRUE)
    diagonals(C) <- sqrt(D)
    s <- slpmvnorm(matrix(lwr), matrix(upr), mean = matrix(mn), chol = C, w = w)
    chk(lpmvnorm(matrix(lwr), matrix(upr), mean = matrix(mn), chol = C, w = w),
        sum(ref$ll), tolerance = 1e-10)
    chk(s$logLik, sum(ref$ll), tolerance = 1e-10)
    for (nm in c("lower", "upper", "mean"))
        chk(c(s[[nm]]), ref[[nm]], tolerance = 1e-8)
}

### J = 1 in the tail against the exact univariate value
lwr <- matrix(c(25, -32, 30, -Inf, -1), nrow = 1)
upr <- matrix(c(27, -30, Inf, -28, 1), nrow = 1)
mn <- matrix(c(.3, -.2, .1, 0, 0), nrow = 1)
for (sd in c(1, 1.7)) {
    C1 <- ltMatrices(sd, diag = TRUE)
    s <- slpmvnorm(lwr, upr, mean = mn, chol = C1)
    ref <- uni(c(lwr), c(upr), c(mn), sd)
    chk(c(lpmvnorm(lwr, upr, mean = mn, chol = C1, logLik = FALSE)), ref$ll,
        tolerance = 1e-10)
    chk(c(s$logLik), ref$ll, tolerance = 1e-10)
    for (nm in c("lower", "upper", "mean"))
        chk(c(s[[nm]]), ref[[nm]], tolerance = 1e-10)
    ### score wrt chol = sd: -(a * dlower + b * dupper) for standardised
    ### finite limits a, b
    a1 <- (c(lwr) - c(mn)) / sd
    b1 <- (c(upr) - c(mn)) / sd
    chk(c(s$chol), -ifelse(is.finite(a1), a1 * ref$lower, 0) -
                    ifelse(is.finite(b1), b1 * ref$upper, 0), tolerance = 1e-10)
}

### 2. scores vs numerical derivatives far in the tails, correlated case

if (require("numDeriv", quietly = TRUE)) {

    ### boxes shifted into the lower / upper / both tails; one chol per
    ### observation, so the gradient of the sum is the per-observation score
    J <- 3; N <- 2
    w <- matrix(runif((J - 1) * 30), nrow = J - 1)
    for (shift in c(10, 35))
    for (sg in list(rep(-1, J), rep(1, J), c(1, -1, 1))) {
        ctr <- sg * shift + matrix(rnorm(J * N, sd = .3), nrow = J)
        wd <- matrix(runif(J * N, .5, 2), nrow = J)
        lwr <- ctr - wd
        upr <- ctr + wd
        mn <- matrix(rnorm(J * N, sd = .3), nrow = J)
        C <- ltMatrices(matrix(runif(J * (J + 1) / 2 * N, -.6, .6), ncol = N),
                        diag = TRUE)
        diagonals(C) <- matrix(runif(J * N, .7, 1.5), ncol = N)
        ll <- function(lower = lwr, upper = upr, mean = mn, prm = unclass(C))
            lpmvnorm(lower, upper, mean = mean, w = w,
                     chol = ltMatrices(matrix(prm, ncol = N), diag = TRUE))
        s <- slpmvnorm(lwr, upr, mean = mn, chol = C, w = w)
        ### below the floor of the old lpmvnorm()
        stopifnot(allfinite(s),
                  all(s$logLik < log(.Machine$double.eps) - log(ncol(w))))
        chk(grad(function(x) ll(lower = x), lwr), c(s$lower), tolerance = 1e-5)
        chk(grad(function(x) ll(upper = x), upr), c(s$upper), tolerance = 1e-5)
        chk(grad(function(x) ll(mean = x), mn), c(s$mean), tolerance = 1e-5)
        chk(grad(function(x) ll(prm = x), unclass(C)), c(unclass(s$chol)),
            tolerance = 1e-5)
    }
}

### 3. degenerate situations

J <- 4; N <- 3
C <- ltMatrices(matrix(runif(J * (J + 1) / 2 * N, -.6, .6), ncol = N),
                diag = TRUE)
diagonals(C) <- matrix(runif(J * N, .7, 1.5), ncol = N)
W <- matrix(runif((J - 1) * 50), nrow = J - 1)

### evaluate both functions, check consistency and finiteness;
### infinite limits have score 0
robust <- function(lwr, upr, w = W) {
    ll <- lpmvnorm(lwr, upr, chol = C, w = w, logLik = FALSE)
    s <- slpmvnorm(lwr, upr, chol = C, w = w)
    chk(s$logLik, ll, tolerance = 1e-12)
    stopifnot(is.finite(ll), allfinite(s),
              all(s$lower[!is.finite(lwr)] == 0),
              all(s$upper[!is.finite(upr)] == 0))
}

### nodes containing exact zeros and ones (as for an unscrambled Sobol
### sequence), with infinite limits, also far in the tails
w0 <- W
w0[, 1] <- 0
w0[, 2] <- 1
w0[1, 3:4] <- 0
lwr <- cbind(c(-Inf, -.5, -2, -Inf), rep(-1, J), rep(-Inf, J))
upr <- cbind(c(Inf, 1, Inf, .3), rep(Inf, J), rep(Inf, J))
for (sh in c(-30, 0, 30))
    robust(lwr + sh, upr + sh, w = w0)

### tol and fast are ignored
lwr <- matrix(25 - runif(J * N), nrow = J)
upr <- lwr + 2
stopifnot(identical(lpmvnorm(lwr, upr, chol = C, w = W),
                    lpmvnorm(lwr, upr, chol = C, w = W, tol = 1e-100, fast = TRUE)),
          identical(slpmvnorm(lwr, upr, chol = C, w = W),
                    slpmvnorm(lwr, upr, chol = C, w = W, tol = 1e-100, fast = TRUE)))

### empty box in observation 2 (lower == upper, finite or +-Inf, in the
### first or last dimension): exactly -Inf and zero scores for this
### observation, the other observations are unaffected
mn <- matrix(rnorm(J * N, sd = .3), nrow = J)
for (d in c(1, J))
for (v in c(.5, -Inf, Inf)) {
    lwr <- -matrix(runif(J * N, .5, 2), nrow = J)
    upr <- matrix(runif(J * N, .5, 2), nrow = J)
    lwr[d, 2] <- upr[d, 2] <- v
    ll <- lpmvnorm(lwr, upr, mean = mn, chol = C, w = W, logLik = FALSE)
    s <- slpmvnorm(lwr, upr, mean = mn, chol = C, w = W)
    chk(s$logLik, ll)
    stopifnot(ll[2] == -Inf, is.finite(ll[-2]), allfinite(s[-1]),
              all(sapply(s[-1], function(x) all(unclass(x)[, 2] == 0))))
}
