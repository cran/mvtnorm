### R code from vignette source 'lmvnorm_src.Rnw'

###################################################
### code chunk number 1: mvtnorm-citation
###################################################
year <- substr(packageDescription("mvtnorm")$Date, 1, 4)
version <- packageDescription("mvtnorm")$Version


###################################################
### code chunk number 2: chk
###################################################
chk <- function(...) stopifnot(isTRUE(all.equal(...)))


###################################################
### code chunk number 3: example
###################################################
library("mvtnorm")
set.seed(290875)
N <- 4L
J <- 5L
rn <- paste0("C_", 1:N)
nm <- LETTERS[1:J]
Jn <- J * (J - 1) / 2
## data
xn <- matrix(runif(N * Jn), ncol = N)
colnames(xn) <- rn
xd <- matrix(runif(N * (Jn + J)), ncol = N)
colnames(xd) <- rn

(lxn <- ltMatrices(xn, byrow = TRUE, names = nm))
dim(lxn)
dimnames(lxn)
lxd <- ltMatrices(xd, byrow = TRUE, diag = TRUE, names = nm)
dim(lxd)
dimnames(lxd)

class(lxn) <- "syMatrices"
lxn


###################################################
### code chunk number 4: ex-reorder
###################################################
## constructor + .reorder + as.array
a <- as.array(ltMatrices(xn, byrow = TRUE))
b <- as.array(ltMatrices(ltMatrices(xn, byrow = TRUE), 
                         byrow = FALSE))
chk(a, b)

a <- as.array(ltMatrices(xn, byrow = FALSE))
b <- as.array(ltMatrices(ltMatrices(xn, byrow = FALSE), 
                         byrow = TRUE))
chk(a, b)

a <- as.array(ltMatrices(xd, byrow = TRUE, diag = TRUE))
b <- as.array(ltMatrices(ltMatrices(xd, byrow = TRUE, diag = TRUE), 
                         byrow = FALSE))
chk(a, b)

a <- as.array(ltMatrices(xd, byrow = FALSE, diag = TRUE))
b <- as.array(ltMatrices(ltMatrices(xd, byrow = FALSE, diag = TRUE), 
                         byrow = TRUE))
chk(a, b)


###################################################
### code chunk number 5: ex-subset
###################################################
## subset
a <- as.array(ltMatrices(xn, byrow = FALSE)[1:2, 2:4])
b <- as.array(ltMatrices(xn, byrow = FALSE))[2:4, 2:4, 1:2]
chk(a, b)

a <- as.array(ltMatrices(xn, byrow = TRUE)[1:2, 2:4])
b <- as.array(ltMatrices(xn, byrow = TRUE))[2:4, 2:4, 1:2]
chk(a, b)

a <- as.array(ltMatrices(xd, byrow = FALSE, diag = TRUE)[1:2, 2:4])
b <- as.array(ltMatrices(xd, byrow = FALSE, diag = TRUE))[2:4, 2:4, 1:2]
chk(a, b)

a <- as.array(ltMatrices(xd, byrow = TRUE, diag = TRUE)[1:2, 2:4])
b <- as.array(ltMatrices(xd, byrow = TRUE, diag = TRUE))[2:4, 2:4, 1:2]
chk(a, b)


###################################################
### code chunk number 6: ex-subset-2
###################################################
## subset
j <- c(1, 3, 5)
a <- as.array(ltMatrices(xn, byrow = FALSE)[1:2, j])
b <- as.array(ltMatrices(xn, byrow = FALSE))[j, j, 1:2]
chk(a, b)

a <- as.array(ltMatrices(xn, byrow = TRUE)[1:2, j])
b <- as.array(ltMatrices(xn, byrow = TRUE))[j, j, 1:2]
chk(a, b)

a <- as.array(ltMatrices(xd, byrow = FALSE, diag = TRUE)[1:2, j])
b <- as.array(ltMatrices(xd, byrow = FALSE, diag = TRUE))[j, j, 1:2]
chk(a, b)

a <- as.array(ltMatrices(xd, byrow = TRUE, diag = TRUE)[1:2, j])
b <- as.array(ltMatrices(xd, byrow = TRUE, diag = TRUE))[j, j, 1:2]
chk(a, b)


###################################################
### code chunk number 7: ex-subset-3
###################################################
## subset
j <- -c(1, 3, 5)
a <- as.array(ltMatrices(xn, byrow = FALSE)[1:2, j])
b <- as.array(ltMatrices(xn, byrow = FALSE))[j, j, 1:2]
chk(a, b)

a <- as.array(ltMatrices(xn, byrow = TRUE)[1:2, j])
b <- as.array(ltMatrices(xn, byrow = TRUE))[j, j, 1:2]
chk(a, b)

a <- as.array(ltMatrices(xd, byrow = FALSE, diag = TRUE)[1:2, j])
b <- as.array(ltMatrices(xd, byrow = FALSE, diag = TRUE))[j, j, 1:2]
chk(a, b)

a <- as.array(ltMatrices(xd, byrow = TRUE, diag = TRUE)[1:2, j])
b <- as.array(ltMatrices(xd, byrow = TRUE, diag = TRUE))[j, j, 1:2]
chk(a, b)


###################################################
### code chunk number 8: ex-subset-4
###################################################
## subset
j <- sample(1:J)
ltM <- ltMatrices(xn, byrow = FALSE)
try(ltM[1:2, j])
class(ltM) <- "syMatrices"
a <- as.array(ltM[1:2, j])
b <- as.array(ltM)[j, j, 1:2]
chk(a, b)


###################################################
### code chunk number 9: ex-Lower_tri
###################################################
## J <- 4
M <- ltMatrices(matrix(1:10, nrow = 10, ncol = 2), diag = TRUE)
Lower_tri(M, diag = FALSE)
Lower_tri(M, diag = TRUE)
M <- ltMatrices(matrix(1:6, nrow = 6, ncol = 2), diag = FALSE)
Lower_tri(M, diag = FALSE)
Lower_tri(M, diag = TRUE)
## multiple symmetric matrices
Lower_tri(invchol2cor(M))


###################################################
### code chunk number 10: ex-diag
###################################################
all(diagonals(ltMatrices(xn, byrow = TRUE)) == 1L)


###################################################
### code chunk number 11: ex-addiag
###################################################
lxd2 <- lxn
diagonals(lxd2) <- 1
chk(as.array(lxd2), as.array(lxn))


###################################################
### code chunk number 12: ex-diagJ
###################################################
(I5 <- diagonals(5L))
diagonals(I5) <- 1:5
I5


###################################################
### code chunk number 13: ex-mult
###################################################
lxn <- ltMatrices(xn, byrow = TRUE)
lxd <- ltMatrices(xd, byrow = TRUE, diag = TRUE)
y <- matrix(runif(N * J), nrow = J)
a <- Mult(lxn, y)
A <- as.array(lxn)
b <- do.call("rbind", lapply(1:ncol(y), 
    function(i) t(A[,,i] %*% y[,i,drop = FALSE])))
chk(a, t(b), check.attributes = FALSE)

a <- Mult(lxd, y)
A <- as.array(lxd)
b <- do.call("rbind", lapply(1:ncol(y), 
    function(i) t(A[,,i] %*% y[,i,drop = FALSE])))
chk(a, t(b), check.attributes = FALSE)

### recycle C
chk(Mult(lxn[rep(1, N),], y), Mult(lxn[1,], y), check.attributes = FALSE)

### recycle y
chk(Mult(lxn, y[,1]), Mult(lxn, y[,rep(1, N)]))

### tcrossprod as multiplication
i <- sample(1:N)[1]
M <- t(as.array(lxn)[,,i])
a <- sapply(1:J, function(j) Mult(lxn[i,], M[,j,drop = FALSE]))
rownames(a) <- colnames(a) <- dimnames(lxn)[[2L]]
b <- as.array(Tcrossprod(lxn[i,]))[,,1]
chk(a, b, check.attributes = FALSE)


###################################################
### code chunk number 14: ex-tmult
###################################################
a <- Mult(lxn, y, transpose = TRUE)
A <- as.array(lxn)
b <- do.call("rbind", lapply(1:ncol(y), 
    function(i) t(t(A[,,i]) %*% y[,i,drop = FALSE])))
chk(a, t(b), check.attributes = FALSE)

a <- Mult(lxd, y, transpose = TRUE)
A <- as.array(lxd)
b <- do.call("rbind", lapply(1:ncol(y), 
    function(i) t(t(A[,,i]) %*% y[,i,drop = FALSE])))
chk(a, t(b), check.attributes = FALSE)

### recycle C
chk(Mult(lxn[rep(1, N),], y, transpose = TRUE), 
    Mult(lxn[1,], y, transpose = TRUE), check.attributes = FALSE)

### recycle y
chk(Mult(lxn, y[,1], transpose = TRUE), 
    Mult(lxn, y[,rep(1, N)], transpose = TRUE))


###################################################
### code chunk number 15: ex-solve
###################################################
## solve
A <- as.array(lxn)
a <- solve(lxn)
a <- as.array(a)
b <- array(apply(A, 3L, function(x) solve(x), simplify = TRUE), 
           dim = rev(dim(lxn)))
chk(a, b, check.attributes = FALSE)

A <- as.array(lxd)
a <- as.array(solve(lxd))
b <- array(apply(A, 3L, function(x) solve(x), simplify = TRUE), 
           dim = rev(dim(lxd)))
chk(a, b, check.attributes = FALSE)

chk(solve(lxn, y), Mult(solve(lxn), y))
chk(solve(lxd, y), Mult(solve(lxd), y))

### recycle C
chk(solve(lxn[1,], y), as.array(solve(lxn[1,]))[,,1] %*% y)
chk(solve(lxn[rep(1, N),], y), solve(lxn[1,], y), check.attributes = FALSE)

### recycle y
chk(solve(lxn, y[,1]), solve(lxn, y[,rep(1, N)]))


###################################################
### code chunk number 16: ex-tsolve
###################################################
chk(solve(lxn[1,], y, transpose = TRUE), 
    t(as.array(solve(lxn[1,]))[,,1]) %*% y)


###################################################
### code chunk number 17: ex-tcrossprod
###################################################
## Tcrossprod
a <- as.array(Tcrossprod(lxn))
b <- array(apply(as.array(lxn), 3L, function(x) tcrossprod(x), simplify = TRUE), 
           dim = rev(dim(lxn)))
chk(a, b, check.attributes = FALSE)

# diagonal elements only
d <- Tcrossprod(lxn, diag_only = TRUE)
chk(d, apply(a, 3, diag))
chk(d, diagonals(Tcrossprod(lxn)))

a <- as.array(Tcrossprod(lxd))
b <- array(apply(as.array(lxd), 3L, function(x) tcrossprod(x), simplify = TRUE), 
           dim = rev(dim(lxd)))
chk(a, b, check.attributes = FALSE)

# diagonal elements only
d <- Tcrossprod(lxd, diag_only = TRUE)
chk(d, apply(a, 3, diag))
chk(d, diagonals(Tcrossprod(lxd)))


###################################################
### code chunk number 18: ex-crossprod
###################################################
## Crossprod
a <- as.array(Crossprod(lxn))
b <- array(apply(as.array(lxn), 3L, function(x) crossprod(x), simplify = TRUE), 
           dim = rev(dim(lxn)))
chk(a, b, check.attributes = FALSE)

# diagonal elements only
d <- Crossprod(lxn, diag_only = TRUE)
chk(d, apply(a, 3, diag))
chk(d, diagonals(Crossprod(lxn)))

a <- as.array(Crossprod(lxd))
b <- array(apply(as.array(lxd), 3L, function(x) crossprod(x), simplify = TRUE), 
           dim = rev(dim(lxd)))
chk(a, b, check.attributes = FALSE)

# diagonal elements only
d <- Crossprod(lxd, diag_only = TRUE)
chk(d, apply(a, 3, diag))
chk(d, diagonals(Crossprod(lxd)))


###################################################
### code chunk number 19: chol
###################################################
Sigma <- Tcrossprod(lxd)
chk(chol(Sigma), lxd)
Sigma <- Tcrossprod(lxn)
## Sigma and chol(Sigma) always have diagonal, lxn doesn't
chk(as.array(chol(Sigma)), as.array(lxn))


###################################################
### code chunk number 20: kronecker
###################################################
J <- 10

d <- TRUE
L <- diag(J)
L[lower.tri(L, diag = d)] <- prm <- runif(J * (J + c(-1, 1)[d + 1]) / 2)

C <- solve(L)

D <- -kronecker(t(C), C)

S <- diag(J)
S[lower.tri(S, diag = TRUE)] <- x <- runif(J * (J + 1) / 2)

SD0 <- matrix(c(S) %*% D, ncol = J)

SD1 <- -crossprod(C, tcrossprod(S, C))

a <- ltMatrices(C[lower.tri(C, diag = TRUE)], diag = TRUE, byrow = FALSE)
b <- ltMatrices(x, diag = TRUE, byrow = FALSE)

SD2 <- -vectrick(a, b, a)
SD2a <- -vectrick(a, b)
chk(SD2, SD2a)

chk(SD0[lower.tri(SD0, diag = d)], 
    SD1[lower.tri(SD1, diag = d)])
chk(SD0[lower.tri(SD0, diag = d)],
    c(unclass(SD2)))

### same; but SD2 is vec(SD0)
S <- t(matrix(as.array(b), byrow = FALSE, nrow = 1))
SD2 <- -vectrick(a, S, a)
SD2a <- -vectrick(a, S)
chk(SD2, SD2a)

chk(c(SD0), c(SD2))

### N > 1
N <- 4L
prm <- runif(J * (J - 1) / 2)
C <- ltMatrices(prm)
S <- matrix(runif(J^2 * N), ncol = N)
A <- vectrick(C, S, C)
Cx <- as.array(C)[,,1]
B <- apply(S, 2, function(x) t(Cx) %*% matrix(x, ncol = J) %*% t(Cx))
chk(A, B)

A <- vectrick(C, S, C, transpose = c(FALSE, FALSE))
Cx <- as.array(C)[,,1]
B <- apply(S, 2, function(x) Cx %*% matrix(x, ncol = J) %*% Cx)
chk(A, B)


###################################################
### code chunk number 21: conv-ex-1
###################################################
prec2pc <- function(x) {
    ret <- -cov2cor(x)
    diag(ret) <- 0
    ret
}
L <- lxn
Sigma <- apply(as.array(L), 3, 
               function(x) tcrossprod(solve(x)), simplify = FALSE)
Prec <- lapply(Sigma, solve)
Corr <- lapply(Sigma, cov2cor)
CP <- lapply(Corr, solve)
PC <- lapply(Prec, function(x) prec2pc(x))
chk(unlist(Sigma), c(as.array(invchol2cov(L))), 
    check.attributes = FALSE)
chk(unlist(Prec), c(as.array(invchol2pre(L))), 
    check.attributes = FALSE)
chk(unlist(Corr), c(as.array(invchol2cor(L))), 
    check.attributes = FALSE)
chk(unlist(CP), c(as.array(Crossprod(invcholD(L)))), 
    check.attributes = FALSE)
chk(unlist(PC), c(as.array(invchol2pc(L))), 
    check.attributes = FALSE)


###################################################
### code chunk number 22: conv-ex-2
###################################################
C <- lxn
Sigma <- apply(as.array(C), 3, 
               function(x) tcrossprod(x), simplify = FALSE)
Prec <- lapply(Sigma, solve)
Corr <- lapply(Sigma, cov2cor)
CP <- lapply(Corr, solve)
PC <- lapply(Prec, function(x) prec2pc(x))
chk(unlist(Sigma), c(as.array(chol2cov(C))), 
    check.attributes = FALSE)
chk(unlist(Prec), c(as.array(chol2pre(C))), 
    check.attributes = FALSE)
chk(unlist(Corr), c(as.array(chol2cor(C))), 
    check.attributes = FALSE)
chk(unlist(CP), c(as.array(Crossprod(solve(Dchol(C))))), 
    check.attributes = FALSE)
chk(unlist(PC), c(as.array(chol2pc(C))), 
    check.attributes = FALSE)


###################################################
### code chunk number 23: conv-ex-3
###################################################
L <- lxd
Sigma <- apply(as.array(L), 3, 
               function(x) tcrossprod(solve(x)), simplify = FALSE)
Prec <- lapply(Sigma, solve)
Corr <- lapply(Sigma, cov2cor)
CP <- lapply(Corr, solve)
PC <- lapply(Prec, function(x) prec2pc(x))
chk(unlist(Sigma), c(as.array(invchol2cov(L))), 
    check.attributes = FALSE)
chk(unlist(Prec), c(as.array(invchol2pre(L))), 
    check.attributes = FALSE)
chk(unlist(Corr), c(as.array(invchol2cor(L))), 
    check.attributes = FALSE)
chk(unlist(CP), c(as.array(Crossprod(invcholD(L)))), 
    check.attributes = FALSE)
chk(unlist(PC), c(as.array(invchol2pc(L))), 
    check.attributes = FALSE)


###################################################
### code chunk number 24: conv-ex-4
###################################################
C <- lxd
Sigma <- apply(as.array(C), 3, 
               function(x) tcrossprod(x), simplify = FALSE)
Prec <- lapply(Sigma, solve)
Corr <- lapply(Sigma, cov2cor)
CP <- lapply(Corr, solve)
PC <- lapply(Prec, function(x) prec2pc(x))
chk(unlist(Sigma), c(as.array(chol2cov(C))), 
    check.attributes = FALSE)
chk(unlist(Prec), c(as.array(chol2pre(C))), 
    check.attributes = FALSE)
chk(unlist(Corr), c(as.array(chol2cor(C))), 
    check.attributes = FALSE)
chk(unlist(CP), c(as.array(Crossprod(solve(Dchol(C))))), 
    check.attributes = FALSE)
chk(unlist(PC), c(as.array(chol2pc(C))), 
    check.attributes = FALSE)


###################################################
### code chunk number 25: aperm-tests
###################################################
L <- lxn
J <- dim(L)[2L]
Lp <- aperm(a = L, perm = p <- sample(1:J), is_chol = FALSE)
chk(invchol2cov(L)[,p], invchol2cov(Lp))

C <- lxn
J <- dim(C)[2L]
Cp <- aperm(a = C, perm = p <- sample(1:J), is_chol = TRUE)
chk(chol2cov(C)[,p], chol2cov(Cp))


###################################################
### code chunk number 26: marg
###################################################
Sigma <- Tcrossprod(lxd)
j <- 1:3
chk(Sigma[,j], Tcrossprod(marg_mvnorm(chol = lxd, which = j)$chol))
j <- 2:4
chk(Sigma[,j], Tcrossprod(marg_mvnorm(chol = lxd, which = j)$chol))

Sigma <- Tcrossprod(solve(lxd))
j <- 1:3
chk(Sigma[,j], Tcrossprod(solve(marg_mvnorm(invchol = lxd, which = j)$invchol)))
j <- 2:4
chk(Sigma[,j], Tcrossprod(solve(marg_mvnorm(invchol = lxd, which = j)$invchol)))


###################################################
### code chunk number 27: cond-general
###################################################
Sigma <- as.array(Tcrossprod(lxd))[,,1]
j <- 2:4
y <- matrix(c(-1, 2, 1), nrow = 3)

cm <- Sigma[-j, j,drop = FALSE] %*% solve(Sigma[j,j]) %*%  y
cS <- Sigma[-j, -j] - Sigma[-j,j,drop = FALSE] %*% 
      solve(Sigma[j,j]) %*% Sigma[j,-j,drop = FALSE]

cmv <- cond_mvnorm(chol = lxd[1,], which = j, given = y)

chk(cm, cmv$mean)
chk(cS, as.array(Tcrossprod(cmv$chol))[,,1])

Sigma <- as.array(Tcrossprod(solve(lxd)))[,,1]
j <- 2:4
y <- matrix(c(-1, 2, 1), nrow = 3)

cm <- Sigma[-j, j,drop = FALSE] %*% solve(Sigma[j,j]) %*%  y
cS <- Sigma[-j, -j] - Sigma[-j,j,drop = FALSE] %*% 
      solve(Sigma[j,j]) %*% Sigma[j,-j,drop = FALSE]

cmv <- cond_mvnorm(invchol = lxd[1,], which = j, given = y)

chk(cm, cmv$mean)
chk(cS, as.array(Tcrossprod(solve(cmv$invchol)))[,,1])


###################################################
### code chunk number 28: cond-simple
###################################################
Sigma <- as.array(Tcrossprod(lxd))[,,1]
j <- 1:3
y <- matrix(c(-1, 2, 1), nrow = 3)

cm <- Sigma[-j, j,drop = FALSE] %*% solve(Sigma[j,j]) %*%  y
cS <- Sigma[-j, -j] - Sigma[-j,j,drop = FALSE] %*% 
      solve(Sigma[j,j]) %*% Sigma[j,-j,drop = FALSE]

cmv <- cond_mvnorm(chol = lxd[1,], which = j, given = y)

chk(c(cm), c(cmv$mean))
chk(cS, as.array(Tcrossprod(cmv$chol))[,,1])

Sigma <- as.array(Tcrossprod(solve(lxd)))[,,1]
j <- 1:3
y <- matrix(c(-1, 2, 1), nrow = 3)

cm <- Sigma[-j, j,drop = FALSE] %*% solve(Sigma[j,j]) %*%  y
cS <- Sigma[-j, -j] - Sigma[-j,j,drop = FALSE] %*% 
      solve(Sigma[j,j]) %*% Sigma[j,-j,drop = FALSE]

cmv <- cond_mvnorm(invchol = lxd[1,], which = j, given = y)

chk(c(cm), c(cmv$mean))
chk(cS, as.array(Tcrossprod(solve(cmv$invchol)))[,,1])


###################################################
### code chunk number 29: ex-MV
###################################################
N <- 1000L
J <- 50L
lt <- ltMatrices(matrix(runif(N * J * (J + 1) / 2) + 1, ncol = N), 
                 diag = TRUE, byrow = FALSE)
Z <- matrix(rnorm(N * J), ncol = N)
Y <- solve(lt, Z)
ll1 <- sum(dnorm(Mult(lt, Y), log = TRUE)) + sum(log(diagonals(lt)))

S <- as.array(Tcrossprod(solve(lt)))
ll2 <- sum(sapply(1:N, function(i) dmvnorm(x = Y[,i], sigma = S[,,i], log = TRUE)))
chk(ll1, ll2)


###################################################
### code chunk number 30: ex-MV-d
###################################################
ll3 <- ldmvnorm(obs = Y, invchol = lt)
chk(ll1, ll3)


###################################################
### code chunk number 31: ex-MV-mc
###################################################
## marginal of and conditional on these
(j <- 1:5 * 10)
md <- marg_mvnorm(invchol = lt, which = j)
cd <- cond_mvnorm(invchol = lt, which = j, given = Y[j,])

ll3 <- sum(dnorm(Mult(md$invchol, Y[j,]), log = TRUE)) + 
       sum(log(diagonals(md$invchol))) +
       sum(dnorm(Mult(cd$invchol, Y[-j,] - cd$mean), log = TRUE)) + 
       sum(log(diagonals(cd$invchol)))
chk(ll1, ll3)


###################################################
### code chunk number 32: chapterseed
###################################################
set.seed(270312)


###################################################
### code chunk number 33: fct-lpmvnormR
###################################################

lpmvnormR <- function(lower, upper, mean = 0, center = NULL, chol, logLik = TRUE, ...) {

    
    if (!is.matrix(lower)) lower <- matrix(lower, ncol = 1)
    if (!is.matrix(upper)) upper <- matrix(upper, ncol = 1)
    stopifnot(isTRUE(all.equal(dim(lower), dim(upper))))

    stopifnot(inherits(chol, "ltMatrices"))
    byrow_orig <- attr(chol, "byrow")
    chol <- ltMatrices(chol, byrow = TRUE)
    d <- dim(chol)
    ### allow single matrix C
    N <- ifelse(d[1L] == 1, ncol(lower), d[1L])
    J <- d[2L]

    stopifnot(nrow(lower) == J && ncol(lower) == N)
    stopifnot(nrow(upper) == J && ncol(upper) == N)
    if (is.matrix(mean))
        stopifnot(nrow(mean) == J && ncol(mean) == N)

    lower <- lower - mean
    upper <- upper - mean

    if (!is.null(center)) {
        if (!is.matrix(center)) center <- matrix(center, ncol = 1)
        stopifnot(nrow(center) == J && ncol(center == N))
    }
    

    sigma <- Tcrossprod(chol)
    S <- as.array(sigma)
    idx <- 1

    ret <- error <- numeric(N)
    for (i in 1:N) {
        if (dim(sigma)[[1L]] > 1) idx <- i
        tmp <- pmvnorm(lower = lower[,i], upper = upper[,i], sigma = S[,,idx], ...)
        ret[i] <- tmp
        error[i] <- attr(tmp, "error")
    }
    attr(ret, "error") <- error

    if (logLik)
        return(sum(log(pmax(ret, .Machine$double.eps))))

    ret
}



###################################################
### code chunk number 34: ex-lpmvnorm_R
###################################################
J <- 5L
N <- 10L

x <- matrix(runif(N * J * (J + 1) / 2), ncol = N)
lx <- ltMatrices(x, byrow = TRUE, diag = TRUE)

a <- matrix(runif(N * J), nrow = J) - 2
a[sample(J * N)[1:2]] <- -Inf
b <- a + 2 + matrix(runif(N * J), nrow = J)
b[sample(J * N)[1:2]] <- Inf

(phat <- c(lpmvnormR(a, b, chol = lx, logLik = FALSE)))


###################################################
### code chunk number 35: ex-again
###################################################
phat
exp(lpmvnorm(a, b, chol = lx, M = 25000, logLik = FALSE, fast = TRUE))
exp(lpmvnorm(a, b, chol = lx, M = 25000, logLik = FALSE, fast = FALSE))


###################################################
### code chunk number 36: ex-lpmvnorm
###################################################
M <- 10000L
if (require("qrng", quietly = TRUE)) {
    ### quasi-Monte-Carlo
    W <- t(ghalton(M, d = J - 1))
} else {
    ### Monte-Carlo
    W <- matrix(runif(M * (J - 1)), nrow = J - 1)
}

### Genz & Bretz, 2001, without early stopping (really?)
pGB <- lpmvnormR(a, b, chol = lx, logLik = FALSE, 
                algorithm = GenzBretz(maxpts = M, abseps = 0, releps = 0))
### Genz 1992 with quasi-Monte-Carlo, fast pnorm
pGqf <- exp(lpmvnorm(a, b, chol = lx, w = W, M = M, logLik = FALSE, 
                     fast = TRUE))
### Genz 1992, original Monte-Carlo, fast pnorm
pGf <- exp(lpmvnorm(a, b, chol = lx, w = NULL, M = M, logLik = FALSE, 
                    fast = TRUE))
### Genz 1992 with quasi-Monte-Carlo, R::pnorm
pGqs <- exp(lpmvnorm(a, b, chol = lx, w = W, M = M, logLik = FALSE, 
                     fast = FALSE))
### Genz 1992, original Monte-Carlo, R::pnorm
pGs <- exp(lpmvnorm(a, b, chol = lx, w = NULL, M = M, logLik = FALSE, 
                    fast = FALSE))

cbind(pGB, pGqf, pGf, pGqs, pGs)


###################################################
### code chunk number 37: ex-uni
###################################################
### test univariate problem
### call pmvnorm
pGB <- lpmvnormR(a[1,,drop = FALSE], b[1,,drop = FALSE], chol = lx[,1], 
                logLik = FALSE, 
                algorithm = GenzBretz(maxpts = M, abseps = 0, releps = 0))
### call lpmvnorm
pGq <- exp(lpmvnorm(a[1,,drop = FALSE], b[1,,drop = FALSE], chol = lx[,1], 
                   logLik = FALSE))
### ground truth
ptr <- pnorm(b[1,] / c(unclass(lx[,1]))) - pnorm(a[1,] / c(unclass(lx[,1])))

cbind(c(ptr), pGB, pGq)


###################################################
### code chunk number 38: ex-score
###################################################
J <- 5L
N <- 4L

S <- crossprod(matrix(runif(J^2), nrow = J))
prm <- t(chol(S))[lower.tri(S, diag = TRUE)]

### define C
mC <- ltMatrices(matrix(prm, ncol = 1), diag = TRUE)

a <- matrix(runif(N * J), nrow = J) - 2
b <- a + 4
a[2,] <- -Inf
b[3,] <- Inf

M <- 10000L
W <- matrix(runif(M * (J - 1)), ncol = M)

lli <- c(lpmvnorm(a, b, chol = mC, w = W, M = M, logLik = FALSE))

fC <- function(prm) {
    C <- ltMatrices(matrix(prm, ncol = 1), diag = TRUE)
    lpmvnorm(a, b, chol = C, w = W, M = M)
}

sC <- slpmvnorm(a, b, chol = mC, w = W, M = M)

chk(lli, sC$logLik)

if (require("numDeriv", quietly = TRUE))
    chk(grad(fC, unclass(mC)), rowSums(unclass(sC$chol)), check.attributes = FALSE)


###################################################
### code chunk number 39: ex-Lscore
###################################################
mL <- solve(mC)

lliL <- c(lpmvnorm(a, b, invchol = mL, w = W, M = M, logLik = FALSE))

chk(lli, lliL)

fL <- function(prm) {
    L <- ltMatrices(matrix(prm, ncol = 1), diag = TRUE)
    lpmvnorm(a, b, invchol = L, w = W, M = M)
}

sL <- slpmvnorm(a, b, invchol = mL, w = W, M = M)

chk(lliL, sL$logLik)

if (require("numDeriv", quietly = TRUE))
    chk(grad(fL, unclass(mL)), rowSums(unclass(sL$invchol)),
        check.attributes = FALSE)


###################################################
### code chunk number 40: ex-uni-score
###################################################
ptr <- pnorm(b[1,] / c(unclass(mC[,1]))) - pnorm(a[1,] / c(unclass(mC[,1])))
log(ptr)
lpmvnorm(a[1,,drop = FALSE], b[1,,drop = FALSE], chol = mC[,1], logLik = FALSE)
lapply(slpmvnorm(a[1,,drop = FALSE], b[1,,drop = FALSE], chol = mC[,1], logLik =
TRUE), unclass)
sd1 <- c(unclass(mC[,1]))
(dnorm(b[1,] / sd1) * b[1,] - dnorm(a[1,] / sd1) * a[1,]) * (-1) / sd1^2 / ptr


###################################################
### code chunk number 41: chapterseed
###################################################
set.seed(110515)


###################################################
### code chunk number 42: ex-ML-dgp
###################################################
J <- 4
R <- diag(J)
R[1,2] <- R[2,1] <- .25
R[1,3] <- R[3,1] <- .5
R[2,4] <- R[4,2] <- .75
(Sigma <- diag(sqrt(1:J / 2)) %*% R %*% diag(sqrt(1:J / 2)))
(C <- t(chol(Sigma)))


###################################################
### code chunk number 43: ex-ML-C
###################################################
prm <- C[lower.tri(C, diag = TRUE)]
lt <- ltMatrices(matrix(prm, ncol = 1L), 
                 diag = TRUE,    ### has diagonal elements
                 byrow = FALSE)  ### prm is column-major
BYROW <- FALSE   ### later checks
lt <- ltMatrices(lt, 
                 byrow = BYROW)   ### convert to row-major
chk(C, as.array(lt)[,,1], check.attributes = FALSE)
chk(Sigma, as.array(Tcrossprod(lt))[,,1], check.attributes = FALSE)


###################################################
### code chunk number 44: ex-ML-data
###################################################
N <- 100L
Z <- matrix(rnorm(N * J), nrow = J)
Y <- Mult(lt, Z) + (mn <- 1:J)


###################################################
### code chunk number 45: ex-ML-mu-vcov
###################################################
rowMeans(Y)
(Shat <- var(t(Y)) * (N - 1) / N)


###################################################
### code chunk number 46: ex-ML-clogLik
###################################################
Yc <- Y - rowMeans(Y)

ll <- function(parm) {
    C <- ltMatrices(parm, diag = TRUE, byrow = BYROW)
    -ldmvnorm(obs = Yc, chol = C)
}

sc <- function(parm) {
    C <- ltMatrices(parm, diag = TRUE, byrow = BYROW)
    -rowSums(unclass(sldmvnorm(obs = Yc, chol = C)$chol))
}


###################################################
### code chunk number 47: ex-ML-const
###################################################
llim <- rep(-Inf, J * (J + 1) / 2)
llim[which(rownames(unclass(lt)) %in% paste(1:J, 1:J, sep = "."))] <- 1e-4


###################################################
### code chunk number 48: ex-ML-c
###################################################
if (BYROW) {
  cML <- chol(Shat)[upper.tri(Shat, diag = TRUE)]
} else {
  cML <- t(chol(Shat))[lower.tri(Shat, diag = TRUE)]
}
ll(cML)
start <- runif(length(cML))
if (require("numDeriv", quietly = TRUE))
    chk(grad(ll, start), sc(start), check.attributes = FALSE)


###################################################
### code chunk number 49: ex-ML-coptim
###################################################
op <- optim(start, fn = ll, gr = sc, method = "L-BFGS-B", 
            lower = llim, control = list(trace = TRUE))
## ML numerically
ltMatrices(op$par, diag = TRUE, byrow = BYROW)
ll(op$par)
## ML analytically
t(chol(Shat))
ll(cML)
## true C matrix
lt


###################################################
### code chunk number 50: ex-ML-cens
###################################################
prb <- 1:9 / 10
sds <- sqrt(diag(Sigma))
ct <- sapply(1:J, function(j) qnorm(prb, mean = mn[j], sd = sds[j])) 
lwr <- upr <- Y
for (j in 1:J) {
    f <- cut(Y[j,], breaks = c(-Inf, ct[,j], Inf))
    lwr[j,] <- c(-Inf, ct[,j])[f]
    upr[j,] <- c(ct[,j], Inf)[f]
}


###################################################
### code chunk number 51: ex-ML-chk (eval = FALSE)
###################################################
## M <- floor(exp(0:25/10) * 1000)
## lGB <- sapply(M, function(m) {
##     st <- system.time(ret <- 
##         lpmvnormR(lwr, upr, mean = mn, chol = lt, algorithm = 
##                   GenzBretz(maxpts = m, abseps = 0, releps = 0)))
##     return(c(st["user.self"], ll = ret))
## })
## lH <- sapply(M, function(m) {
##     W <- NULL
##     if (require("qrng", quietly = TRUE))
##         W <- t(ghalton(m, d = J - 1))
##     st <- system.time(ret <- lpmvnorm(lwr, upr, mean = mn, 
##                                       chol = lt, w = W, M = m))
##     return(c(st["user.self"], ll = ret))
## })
## lHf <- sapply(M, function(m) {
##     W <- NULL
##     if (require("qrng", quietly = TRUE))
##         W <- t(ghalton(m, d = J - 1))
##     st <- system.time(ret <- lpmvnorm(lwr, upr, mean = mn, chol = lt, 
##                                       w = W, M = m, fast = TRUE))
##     return(c(st["user.self"], ll = ret))
## })


###################################################
### code chunk number 52: ex-ML-fig-data
###################################################
### use pre-computed data, otherwise CRAN complains.
M <-
c(1000, 1105, 1221, 1349, 1491, 1648, 1822, 2013, 2225, 2459, 
2718, 3004, 3320, 3669, 4055, 4481, 4953, 5473, 6049, 6685, 7389, 
8166, 9025, 9974, 11023, 12182)
lGB <-
structure(c(0.054000000000002046, -880.49261248615869, 0.054000000000002046, 
-880.49242591762345, 0.054000000000002046, -880.49299587719224, 
0.053999999999994941, -880.49262902245289, 0.054000000000002046, 
-880.49023141833743, 0.054999999999999716, -880.49278384818467, 
0.054000000000002046, -880.49263211493064, 0.054999999999999716, 
-880.48929666151639, 0.054000000000002046, -880.49251644257947, 
0.054000000000002046, -880.49133879555245, 0.054000000000002046, 
-880.49209059546854, 0.10999999999999943, -880.49160057988672, 
0.11399999999999721, -880.49355264445524, 0.11100000000000421, 
-880.49125049308975, 0.10799999999999699, -880.49215052698423, 
0.10800000000000409, -880.49227539674052, 0.10900000000000176, 
-880.4918790945195, 0.10900000000000176, -880.49200824376248, 
0.19200000000000017, -880.49213191772731, 0.19500000000000028, 
-880.49183887878723, 0.19400000000000261, -880.49213912969094, 
0.19399999999999551, -880.49104161510115, 0.1980000000000004, 
-880.49219786133153, 0.32800000000000296, -880.49160034162105, 
0.3230000000000004, -880.49194070754265, 0.3230000000000004, 
-880.49169841069408), dim = c(2L, 26L), dimnames = list(c("user.self", 
"ll"), NULL))
lH <-
structure(c(0.022999999999996135, -880.48029630630356, 0.027000000000001023, 
-880.49616556706735, 0.028999999999996362, -880.48868341221828, 
0.031999999999996476, -880.49617133610923, 0.034999999999996589, 
-880.48559676487218, 0.039000000000001478, -880.49133269164929, 
0.042999999999999261, -880.49455705455455, 0.047999999999994714, 
-880.49542914386404, 0.052999999999997272, -880.49439060601685, 
0.059999999999995168, -880.48554604736751, 0.064000000000000057, 
-880.49145473103829, 0.071000000000005059, -880.49413777944255, 
0.079000000000000625, -880.49161893244116, 0.087000000000003297, 
-880.49339336832497, 0.094999999999998863, -880.49254091532248, 
0.10599999999999454, -880.49164929759263, 0.117999999999995, 
-880.49250821826354, 0.12900000000000489, -880.49255823649469, 
0.14100000000000534, -880.49250899551407, 0.15699999999999648, 
-880.49044836032715, 0.17300000000000182, -880.4916864321558, 
0.19299999999999784, -880.49117821390132, 0.21099999999999852, 
-880.49228594426586, 0.23299999999999699, -880.49151053053481, 
0.25800000000000267, -880.49153034586084, 0.28699999999999903, 
-880.49192862086477), dim = c(2L, 26L), dimnames = list(c("user.self", 
"ll"), NULL))
lHf <-
structure(c(0.018000000000000682, -880.48706686434321, 0.019000000000005457, 
-880.48863853791863, 0.021999999999998465, -880.48856920618675, 
0.023999999999993804, -880.49392985529778, 0.026000000000003354, 
-880.48602935577401, 0.029000000000003467, -880.49156265644069, 
0.033000000000001251, -880.49941537079019, 0.035000000000003695, 
-880.49445708416158, 0.038000000000003809, -880.49395360381868, 
0.042999999999999261, -880.49364771019248, 0.04700000000000415, 
-880.49295526291712, 0.051999999999999602, -880.4946669985851, 
0.058999999999997499, -880.49374458475813, 0.064999999999997726, 
-880.49419536692346, 0.070000000000000284, -880.49332997377735, 
0.07799999999999585, -880.49145097709174, 0.086000000000005627, 
-880.49237859905554, 0.094000000000001194, -880.49039213464482, 
0.10599999999999454, -880.49106109560728, 0.11500000000000199, 
-880.49157721254892, 0.12899999999999778, -880.49252340367264, 
0.14199999999999591, -880.49102722425221, 0.15800000000000125, 
-880.49208613066503, 0.17099999999999227, -880.4920686898173, 
0.18899999999999295, -880.4922510767417, 0.20799999999999841, 
-880.4923471956962), dim = c(2L, 26L), dimnames = list(c("user.self", 
"ll"), NULL))


###################################################
### code chunk number 53: ex-ML-fig
###################################################
layout(matrix(1:2, nrow = 1))
plot(M, lGB["ll",], ylim = range(c(lGB["ll",], lH["ll",], lHf["ll",])), ylab = "Log-likelihood")
points(M, lH["ll",], pch = 4)
points(M, lHf["ll",], pch = 5)
plot(M, lGB["user.self",], ylim = c(0, max(lGB["user.self",])), ylab = "Time (in sec)")
points(M, lH["user.self",], pch = 4)
points(M, lHf["user.self",], pch = 5)
legend("bottomright", legend = c("pmvnorm", "lpmvnorm", "lpmvnorm(fast)"), pch = c(1, 4, 5), bty = "n")


###################################################
### code chunk number 54: ex-ML-ll
###################################################
M <- 500 
if (require("qrng", quietly = TRUE)) {
    ### quasi-Monte-Carlo
    W <- t(ghalton(M, d = J - 1))
} else {
    ### Monte-Carlo
    W <- matrix(runif(M * (J - 1)), nrow = J - 1)
}
ll <- function(parm, J) {
     m <- parm[1:J]		### mean parameters
     parm <- parm[-(1:J)]	### chol parameters
     C <- matrix(c(parm), ncol = 1L)
     C <- ltMatrices(C, diag = TRUE, byrow = BYROW)
     -lpmvnorm(lower = lwr, upper = upr, mean = m, chol = C, 
               w = W, M = M, logLik = TRUE)
}


###################################################
### code chunk number 55: ex-ML-check
###################################################
prm <- c(mn, unclass(lt))
ll(prm, J = J)
lpmvnormR(lwr, upr, mean = mn, chol = lt, 
         algorithm = GenzBretz(maxpts = M, abseps = 0, releps = 0))
(llprm <- lpmvnorm(lwr, upr, mean = mn, chol = lt, w = W, M = M))
chk(llprm, sum(lpmvnorm(lwr, upr, mean = mn, chol = lt, w = W, 
                        M = M, logLik = FALSE)))


###################################################
### code chunk number 56: ex-ML-sc
###################################################
sc <- function(parm, J) {
    m <- parm[1:J]             ### mean parameters
    parm <- parm[-(1:J)]       ### chol parameters
    C <- matrix(c(parm), ncol = 1L)
    C <- ltMatrices(C, diag = TRUE, byrow = BYROW)
    ret <- slpmvnorm(lower = lwr, upper = upr, mean = m, chol = C, 
                     w = W, M = M, logLik = TRUE)
    return(-c(rowSums(ret$mean), rowSums(unclass(ret$chol))))
}


###################################################
### code chunk number 57: ex-ML-sc-chk
###################################################
if (require("numDeriv", quietly = TRUE))
    chk(grad(ll, prm, J = J), sc(prm, J = J), check.attributes = FALSE)


###################################################
### code chunk number 58: ex-ML
###################################################
llim <- rep(-Inf, J + J * (J + 1) / 2)
llim[J + which(rownames(unclass(lt)) %in% paste(1:J, 1:J, sep = "."))] <- 1e-4

if (BYROW) {
  start <- c(rowMeans(Y), chol(Shat)[upper.tri(Shat, diag = TRUE)])
} else {
  start <- c(rowMeans(Y), t(chol(Shat))[lower.tri(Shat, diag = TRUE)])
}

ll(start, J = J)

op <- optim(start, fn = ll, gr = sc, J = J, method = "L-BFGS-B", 
            lower = llim, control = list(trace = TRUE))

op$value ## compare with 
ll(prm, J = J)


###################################################
### code chunk number 59: ex-ML-C
###################################################
(C <- ltMatrices(matrix(op$par[-(1:J)], ncol = 1), 
                 diag = TRUE, byrow = BYROW))
lt


###################################################
### code chunk number 60: ex-ML-mu
###################################################
op$par[1:J]
mn


###################################################
### code chunk number 61: ex-ML-Shat
###################################################
Tcrossprod(lt)		### true Sigma
Tcrossprod(C)           ### interval-censored obs
Shat                    ### "exact" obs


###################################################
### code chunk number 62: regressions
###################################################
c(cond_mvnorm(chol = C, which = 2:J, given = diag(J - 1))$mean)


###################################################
### code chunk number 63: lm-ex
###################################################
dY <- as.data.frame(t(Y))
colnames(dY) <- paste0("Y", 1:J)
coef(m1 <- lm(Y1 ~ ., data = dY))[-1L]


###################################################
### code chunk number 64: hessian
###################################################
H <- optim(op$par, fn = ll, gr = sc, J = J, method = "L-BFGS-B", 
           lower = llim, hessian = TRUE)$hessian


###################################################
### code chunk number 65: ML-sample
###################################################
L <- t(chol(H))
L <- ltMatrices(L[lower.tri(L, diag = TRUE)], diag = TRUE)
Nsim <- 50000
Z <- matrix(rnorm(Nsim * nrow(H)), ncol = Nsim)
rC <- solve(L, Z)[-(1:J),] + op$par[-(1:J)] ### remove mean parameters


###################################################
### code chunk number 66: ML-check
###################################################
c(sqrt(rowMeans((rC - rowMeans(rC))^2)))
c(sqrt(diagonals(Crossprod(solve(L)))))


###################################################
### code chunk number 67: rC
###################################################
rC <- ltMatrices(rC, diag = TRUE)


###################################################
### code chunk number 68: ML-beta
###################################################
rbeta <- cond_mvnorm(chol = rC, which = 2:J, given = diag(J - 1))$mean
sqrt(rowMeans((rbeta - rowMeans(rbeta))^2))


###################################################
### code chunk number 69: se-ex
###################################################
sqrt(diag(vcov(m1)))[-1L]


###################################################
### code chunk number 70: ex-ML-cd
###################################################
ll_cd <- function(parm, J) {
     m <- parm[1:J]             ### mean parameters
     parm <- parm[-(1:J)]       ### chol parameters
     C <- matrix(c(parm), ncol = 1L)
     C <- ltMatrices(C, diag = TRUE, byrow = BYROW)
     -ldpmvnorm(obs = Y[1:2,], lower = lwr[-(1:2),], 
                upper = upr[-(1:2),], mean = m, chol = C, 
                w = W[-(1:2),,drop = FALSE], M = M)
}
sc_cd <- function(parm, J) {
    m <- parm[1:J]             ### mean parameters
    parm <- parm[-(1:J)]       ### chol parameters
    C <- matrix(c(parm), ncol = 1L)
    C <- ltMatrices(C, diag = TRUE, byrow = BYROW)
    ret <- sldpmvnorm(obs = Y[1:2,], lower = lwr[-(1:2),],
                      upper = upr[-(1:2),], mean = m, chol = C, 
                      w = W[-(1:2),,drop = FALSE], M = M)
    return(-c(rowSums(ret$mean), rowSums(unclass(ret$chol))))
}


###################################################
### code chunk number 71: ex-ML-cd-score
###################################################
if (require("numDeriv", quietly = TRUE))
    chk(grad(ll_cd, start, J = J), sc_cd(start, J = J), 
        check.attributes = FALSE, tol = 1e-6)


###################################################
### code chunk number 72: ex-ML-cd-optim
###################################################
op <- optim(start, fn = ll_cd, gr = sc_cd, J = J, 
            method = "L-BFGS-B", lower = llim, 
            control = list(trace = TRUE))
## estimated C
ltMatrices(matrix(op$par[-(1:J)], ncol = 1), 
           diag = TRUE, byrow = BYROW)
## compare with true C
lt
## estimated means
op$par[1:J]
## compare with true means
mn


###################################################
### code chunk number 73: ex-stand
###################################################
C <- ltMatrices(runif(10))
all.equal(as.array(chol2cov(standardize(chol = C))),
          as.array(chol2cor(standardize(chol = C))))
L <- solve(C)
all.equal(as.array(invchol2cov(standardize(invchol = L))),
          as.array(invchol2cor(standardize(invchol = L))))


###################################################
### code chunk number 74: gc-classical
###################################################
data("iris")
J <- 4
Z <- t(qnorm(do.call("cbind", lapply(iris[1:J], rank)) / (nrow(iris) + 1)))
(CR <- cor(t(Z)))
ll <- function(parm) {
    C <- ltMatrices(parm)
    Cs <- standardize(C)
    -ldmvnorm(obs = Z, chol = Cs)
}
sc <- function(parm) {
    C <- ltMatrices(parm)
    Cs <- standardize(C)
    -rowSums(Lower_tri(destandardize(chol = C, 
        score_schol = sldmvnorm(obs = Z, chol = Cs)$chol)))
}
start <- t(chol(CR))
start <- start[lower.tri(start)]
if (require("numDeriv", quietly = TRUE))
    chk(grad(ll, start), sc(start), check.attributes = FALSE)
op <- optim(start, fn = ll, gr = sc, method = "BFGS", hessian = TRUE)
op$value
S_ML <- chol2cov(standardize(ltMatrices(op$par)))


###################################################
### code chunk number 75: gc-NPML
###################################################
lwr <- do.call("cbind", lapply(iris[1:J], rank, ties.method = "min")) - 1L
upr <- do.call("cbind", lapply(iris[1:J], rank, ties.method = "max"))
lwr <- t(qnorm(lwr / nrow(iris)))
upr <- t(qnorm(upr / nrow(iris)))

M <- 500 
if (require("qrng", quietly = TRUE)) {
    ### quasi-Monte-Carlo
    W <- t(ghalton(M, d = J - 1))
} else {
    ### Monte-Carlo
    W <- matrix(runif(M * (J - 1)), nrow = J - 1)
}

ll <- function(parm) {
    C <- ltMatrices(parm)
    Cs <- standardize(C)
    -lpmvnorm(lower = lwr, upper = upr, chol = Cs, M = M, w = W)
}
sc <- function(parm) {
    C <- ltMatrices(parm)
    Cs <- standardize(C)
    -rowSums(Lower_tri(destandardize(chol = C, 
        score_schol = slpmvnorm(lower = lwr, upper = upr, chol = Cs, 
                               M = M, w = W)$chol)))
}
if (require("numDeriv", quietly = TRUE))
    chk(grad(ll, start), sc(start), check.attributes = FALSE)
op2 <- optim(start, fn = ll, gr = sc, method = "BFGS", hessian = TRUE)
S_NPML <- chol2cov(standardize(ltMatrices(op2$par)))


###################################################
### code chunk number 76: gc
###################################################
S_ML
S_NPML


###################################################
### code chunk number 77: gc-se
###################################################
sd_ML <- ltMatrices(sqrt(diag(solve(op$hessian))))
diagonals(sd_ML) <- 0
sd_NPML <- try(ltMatrices(sqrt(diag(solve(op2$hessian)))))
if (!inherits(sd_NPML, "try-error")) {
    diagonals(sd_NPML) <- 0
    print(sd_ML)
    print(sd_NPML)
}


