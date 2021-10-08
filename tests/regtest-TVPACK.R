library("mvtnorm")

## Showing the TVPACK() gives *NON*-random results:
(cor1 <- toeplitz(c(1, 1/4, -1/8)))
(up1  <- c(1/4, 7/4, 5/8))
d <- length(up1) # = 3
pmvt.. <- function(df, algorithm)
    vapply(df, function(df) pmvt(upper=up1, corr=cor1, df=df, algorithm=algorithm),
           numeric(1))

dfs <- 1:9
pmvt_TV.7 <- replicate(7, pmvt..(dfs, TVPACK()))

stopifnot(pmvt_TV.7 == pmvt_TV.7[,1])
(pmvt.TV. <- pmvt_TV.7[,1])
(pmvt.TV  <- pmvt..(dfs, TVPACK(1e-14)))# has no effect here
all.equal(max(abs(pmvt.TV - pmvt.TV.)), 0) ## all 0 {unexpectedly ??}


set.seed(47) ## and default algorithm: -> *random* result
pmvt_7 <- replicate(7, vapply(dfs, function(df) pmvt(df=df, upper=up1, corr=cor1), numeric(1)))
## relative errors
relE <- 1 - pmvt_7 / pmvt.TV
rng.rE <- range(abs(relE))
stopifnot(1e-6 < rng.rE[1], rng.rE[2] < 7e-4)
stopifnot(all.equal(
    colMeans(abs(relE)),
    c(88, 64, 105, 73, 52, 90, 87)*1e-6, tol= 1e-3))


set.seed(29)

########################################################################
## 3 dim example
corr <- cov2cor(crossprod(matrix(runif(9,-1,1),3,3))+diag(3))
df <- rpois(1,3)+1

## central t distribution (-Inf,upper)
ctrl <- GenzBretz(maxpts = 2500000, abseps = 0.000001, releps = 0)
upper <- rexp(3,1)
pmvt(upper=upper, corr=corr, df = df, algorithm = ctrl)
pmvt(upper=upper, corr=corr, df = df, algorithm = TVPACK())

## central t distribution (lower,Inf)
lower <- -rexp(3,1)
pmvt(lower=lower, upper=rep(Inf,3), corr=corr, df = df, algorithm = ctrl)
pmvt(lower=lower, upper=rep(Inf,3), corr=corr, df = df, algorithm = TVPACK())

## non-central t (not possible for TVPACK)
delt <- rexp(3,1/10)
upper <- delt+runif(3)
ctrl <- GenzBretz(maxpts = 2500000, abseps = 0.000001, releps = 0)
pmvt(upper=upper, corr=corr, df = df, algorithm = ctrl, delta = delt)
tools::assertError(pmvt(upper=upper, corr=corr, df = df, algorithm = TVPACK(), delta = delt))

## central mvn (-Inf, upper)
upper <- rexp(3,1)
pmvnorm(upper=upper, corr=corr, algorithm = ctrl)
pmvnorm(upper=upper, corr=corr, algorithm = TVPACK())

## central mvn (lower, Inf)
lower <- rexp(3,5)
pmvnorm(lower=lower,upper=rep(Inf, 3), corr=corr, algorithm = ctrl)
pmvnorm(lower=lower,upper=rep(Inf, 3), corr=corr, algorithm = TVPACK())

## non-central mvn
delt <- rexp(3,1/10)
upper <- delt+rexp(3,1)
pmvnorm(upper=upper, corr=corr, algorithm = ctrl,     mean = delt)
pmvnorm(upper=upper, corr=corr, algorithm = TVPACK(), mean = delt) # should not error

########################################################################
## 2 dim example
corr <- cov2cor(crossprod(matrix(runif(4,-1,1),2,2))+diag(2))
upper <- rexp(2,1)
df <- rpois(1, runif(1, 0, 20))

## central t (-Inf, upper)
pmvt(upper=upper, corr=corr, df = df, algorithm = ctrl)
pmvt(upper=upper, corr=corr, df = df, algorithm = TVPACK())

## central t (lower, Inf)
pmvt(lower=-upper, upper=rep(Inf, 2), corr=corr, df = df, algorithm = ctrl)
pmvt(lower=-upper, upper=rep(Inf, 2), corr=corr, df = df, algorithm = TVPACK())

## non-central t
delt <- rexp(2,1/5)
upper <- delt+rexp(2,1)
pmvnorm(upper=upper, corr=corr, algorithm = ctrl, mean = delt)
pmvnorm(upper=upper, corr=corr, algorithm = TVPACK(), mean = delt)

########################################################################
## comparison with Miwa
## 2d
corr <- cov2cor(crossprod(matrix(runif(4,-1,1),2,2))+diag(2))
upper <- rexp(2, 1)

pmvnorm(upper=upper, corr=corr, algorithm = Miwa(steps=128))
pmvnorm(upper=upper, corr=corr, algorithm = TVPACK())

## 3d
corr <- cov2cor(crossprod(matrix(runif(9,-1,1),3,3))+diag(3))
upper <- rexp(3, 1)

ctrl <- Miwa(steps=128)
pmvnorm(upper=upper, corr=corr, algorithm = ctrl)
pmvnorm(upper=upper, corr=corr, algorithm = TVPACK())

##==== Cases where some  (lower[j], upper[j]) == (-Inf, Inf) :
S <- toeplitz(c(1, 1/2, 1/4))

set.seed(11)
P0 <- pmvnorm(lower=c(-Inf, 0, 0), upper=Inf, corr=S)
P1 <- pmvnorm(lower=c(-Inf, 0, 0), upper=Inf, corr=S, algorithm = TVPACK()) # had failed
P2 <- pmvnorm(lower=c(-Inf, 0, 0), upper=Inf, corr=S, algorithm = Miwa())
P2a<- pmvnorm(lower=c(-Inf, 0, 0), upper=Inf, corr=S, algorithm = Miwa(512))
P2.<- pmvnorm(lower=c(-Inf, 0, 0), upper=Inf, corr=S, algorithm = Miwa(2048))

stopifnot(all.equal(1/3, c(P0), tol=1e-14)
        , all.equal(1/3, c(P1), tol=1e-14)
        , all.equal(1/3, c(P2), tol=1e-9 ) # 3.765e-10
        , all.equal(1/3, c(P2a),tol=4e-12) # 8.32 e-13
        , all.equal(1/3, c(P2.),tol=2e-12) # 5.28 e-13
)

## t-dist [TVPACK() had failed] :
set.seed(11)
Ptdef <- replicate(20, c(pmvt(lower=c(-Inf, 1, 2), upper=Inf, df=2, corr=S)))
unique(Ptdef)# see length 1; i.e., same result [even though default is Monte-Carlo ??]
Pt1 <- pmvt(lower=c(-Inf, 1, 2), upper=Inf, df=2, corr=S, algorithm = TVPACK())
P. <- 0.0570404044526986
stopifnot(exprs = {
    all.equal(P., c(Pt1), tol = 1e-14)# seen 3.65 e-16
    abs(P. - Ptdef) < 1e-15 # seen 1.39 e-17
})
