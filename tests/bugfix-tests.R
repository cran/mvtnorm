invisible(options(echo = TRUE))
library("mvtnorm")
set.seed(290875)

# correlation matrices for unequal variances were wrong
# from Pamela Ohman-Strickland <ohmanpa@UMDNJ.EDU>

a <- 4.048
shi <- -9
slo <- -10
mu <- -5
sig <- matrix(c(1,1,1,2),ncol=2) 
pmvnorm(lower=c(-a,slo),upper=c(a,shi),mean=c(mu,2*mu),sigma=sig)

# check if set.seed works (starting from 0.5-7)
n <- 5
lower <- -1
upper <- 3
df <- 4
corr <- diag(5)
corr[lower.tri(corr)] <- 0.5
delta <- rep(0, 5)
set.seed(290875)
prob1 <- pmvt(lower=lower, upper=upper, delta=delta, df=df, corr=corr)
set.seed(290875)
prob2 <- pmvt(lower=lower, upper=upper, delta=delta, df=df, corr=corr)
stopifnot(all.equal(prob1, prob2))

# confusion for univariate probabilities when sigma is a matrix
# by Jerome Asselin <jerome@hivnet.ubc.ca>
a <- pmvnorm(lower=-Inf,upper=2,mean=0,sigma=matrix(1.5))
attributes(a) <- NULL
stopifnot(all.equal(a, pnorm(2, sd=sqrt(1.5))))
a <- pmvnorm(lower=-Inf,upper=2,mean=0,sigma=matrix(.5))
attributes(a) <- NULL
stopifnot(all.equal(a, pnorm(2, sd=sqrt(.5))))
a <- pmvnorm(lower=-Inf,upper=2,mean=0,sigma=.5)
attributes(a) <- NULL
stopifnot(all.equal(a, pnorm(2, sd=sqrt(.5))))

# log argument added by Jerome Asselin <jerome@hivnet.ubc.ca>
dmvnorm(x=c(0,0), mean=c(1,1),log=TRUE)
dmvnorm(x=c(0,0), mean=c(25,25),log=TRUE)
dmvnorm(x=c(0,0), mean=c(30,30),log=TRUE)
stopifnot(all.equal(dmvnorm(x=0, mean=30,log=TRUE), 
                    dnorm(0,30,log=TRUE)))

# large df
pnorm(2)^2
pmvt(lower=c(-Inf,-Inf), upper=c(2,2), delta=c(0, 0), df=25, corr=diag(2))
pmvt(lower=c(-Inf,-Inf), upper=c(2,2), delta=c(0, 0), df=250, corr=diag(2))
pmvt(lower=c(-Inf,-Inf), upper=c(2,2), delta=c(0, 0), df=1340, corr=diag(2))
pmvt(lower=c(-Inf,-Inf), upper=c(2,2), delta=c(0, 0), df=2500, corr=diag(2))
pmvt(lower=c(-100,-100), upper=c(2,2), delta=c(0, 0), df=2500, corr=diag(2))

# df = 0
pmvt(lower=c(-Inf,-Inf), upper=c(2,2), delta=c(0, 0), df=0, corr=diag(2))
pmvt(lower=-Inf, upper = 2, delta=0, df=0, corr=1)
pnorm(2)

# larger dimensions
pnorm(2)^2
pmvnorm(lower=rep(-Inf, 2), upper=rep(2,2), sigma = diag(2))
pnorm(2)^90
pmvnorm(lower=rep(-Inf, 90), upper=rep(2,90), sigma = diag(90))
pnorm(2)^199
pmvnorm(lower=rep(-Inf, 199), upper=rep(2,199), sigma = diag(199))

# larger dimensions, again. Spotted by Chihiro Kuroki <kuroki@oak.dti.ne.jp>
# Alan's fix to MVCHNC solves this problem
cr = matrix(0.5, nr = 4, nc = 4)
diag(cr) = 1
cr
a <- pmvt(low = -rep(1, 4), upp = rep(1, 4), df = 999, corr = cr)
b <- pmvt(low = -rep(1, 4), upp = rep(1, 4), df = 4999, corr = cr)
b
attributes(a) <- NULL
attributes(b) <- NULL
stopifnot(all.equal(round(a, 3), round(b, 3)))

# cases where the support is the empty set tried to compute something.
# spotted by Peter Thomson <peter@statsresearch.co.nz>
stopifnot(pmvnorm(upper=c(-Inf,1)) == 0)
stopifnot(pmvnorm(lower=c(Inf,1)) == 0)
stopifnot(pmvnorm(lower=c(-2,0),upper=c(-1,1),corr=matrix(rep(1,4),2,2)) == 0)

# bugged Fritz (long time ago)
stopifnot(all.equal(pmvnorm(-Inf, c(Inf, 0), 0, diag(2)), pmvnorm(-Inf,
                    c(Inf, 0), 0)))

# this is a bug in `mvtdst' nobody was able to fix yet :-(
stopifnot(pmvnorm(lo=c(-Inf,-Inf), up=c(Inf,Inf), mean=c(0,0)) == 1)

### check for correct random seed initialization
### problem reported by Karen Conneely <conneely@umich.edu>
dm <- 250000
iters <- 2
corr <- .7
dim <- 100
abserr <- .0000035
cutoff <- -5.199338
mn <- rep(0,dim)
mat <- diag(dim)
for (i in 1:dim) {
    for (j in 1:(i-1)) {
        mat[i,j]=mat[j,i]=corr^(i-j)
    }
}
ll <- rep(cutoff, dim)
mn <- rep(0, dim)
p <- matrix(0, iters,1)

set.seed(290875)
for (i in 1:iters) {
   pp <- pmvnorm(lower=ll, sigma=mat, maxpts=dm, abseps=abserr)
   p[i] <- 1-pp
}
stopifnot(abs(p[1] - p[2]) < 2 * abserr)
ptmp <- p
set.seed(290875)
for (i in 1:iters) {
   pp <- pmvnorm(lower=ll, sigma=mat, maxpts=dm, abseps=abserr)
   p[i] <- 1-pp
}
stopifnot(all.equal(p, ptmp))

### same for algoritm = Miwa

pmvnormM <- function(...) pmvnorm(..., algorithm = Miwa())

a <- 4.048
shi <- -9
slo <- -10
mu <- -5
sig <- matrix(c(1,1,1,2),ncol=2) 
pmvnormM(lower=c(-a,slo),upper=c(a,shi),mean=c(mu,2*mu),sigma=sig)

# check if set.seed works (starting from 0.5-7)
n <- 5
lower <- -1
upper <- 3
df <- 4
corr <- diag(5)
corr[lower.tri(corr)] <- 0.5
delta <- rep(0, 5)
set.seed(290875)
prob1 <- pmvnormM(lower=lower, upper=upper, mean = delta, corr=corr)
set.seed(290875)
prob2 <- pmvnormM(lower=lower, upper=upper, mean = delta, corr=corr)
stopifnot(all.equal(prob1, prob2))

# confusion for univariate probabilities when sigma is a matrix
# by Jerome Asselin <jerome@hivnet.ubc.ca>
a <- pmvnormM(lower=-Inf,upper=2,mean=0,sigma=matrix(1.5))
attributes(a) <- NULL
stopifnot(all.equal(a, pnorm(2, sd=sqrt(1.5))))
a <- pmvnormM(lower=-Inf,upper=2,mean=0,sigma=matrix(.5))
attributes(a) <- NULL
stopifnot(all.equal(a, pnorm(2, sd=sqrt(.5))))
a <- pmvnormM(lower=-Inf,upper=2,mean=0,sigma=.5)
attributes(a) <- NULL
stopifnot(all.equal(a, pnorm(2, sd=sqrt(.5))))


# cases where the support is the empty set tried to compute something.
# spotted by Peter Thomson <peter@statsresearch.co.nz>
stopifnot(pmvnormM(upper=c(-Inf,1)) == 0)
stopifnot(pmvnormM(lower=c(Inf,1)) == 0)

# bugged Fritz (long time ago)
stopifnot(all.equal(pmvnormM(-Inf, c(Inf, 0), 0, diag(2)), pmvnormM(-Inf,
                    c(Inf, 0), 0)))

# this is a bug in `mvtdst' nobody was able to fix yet :-(
stopifnot(pmvnormM(lo=c(-Inf,-Inf), up=c(Inf,Inf), mean=c(0,0)) == 1)

### check for correct random seed initialization
### problem reported by Karen Conneely <conneely@umich.edu>
dm <- 250000
iters <- 2
corr <- .7
dim <- 10
abserr <- .0000035
cutoff <- -5.199338
mn <- rep(0,dim)
mat <- diag(dim)
for (i in 1:dim) {
    for (j in 1:(i-1)) {
        mat[i,j]=mat[j,i]=corr^(i-j)
    }
}
ll <- rep(cutoff, dim)
mn <- rep(0, dim)
p <- matrix(0, iters,1)

set.seed(290875)
for (i in 1:iters) {
   pp <- pmvnormM(lower=ll, sigma=mat, maxpts=dm, abseps=abserr)
   p[i] <- 1-pp
}
stopifnot(abs(p[1] - p[2]) < 2 * abserr)
ptmp <- p
set.seed(290875)
for (i in 1:iters) {
   pp <- pmvnormM(lower=ll, sigma=mat, maxpts=dm, abseps=abserr)
   p[i] <- 1-pp
}
stopifnot(all.equal(p, ptmp))
