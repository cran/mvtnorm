invisible(options(echo = TRUE))
library(mvtnorm)
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

