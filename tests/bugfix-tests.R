invisible(options(echo = TRUE))
library(mvtnorm)

# correlation matrices for unequal variances were wrong
# from Pamela Ohman-Strickland <ohmanpa@UMDNJ.EDU>

a <- 4.048
shi <- -9
slo <- -10
mu <- -5
sig <- matrix(c(1,1,1,2),ncol=2) 
pmvnorm(lower=c(-a,slo),upper=c(a,shi),mean=c(mu,2*mu),sigma=sig)
