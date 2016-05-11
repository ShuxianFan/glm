
rm(list=ls())

logit <- function(x){
	log(x/(1-x))
}

expit <- function(x){
	exp(x)/(1+exp(x))
}

##########################################################
### Simulate data for binomial GLM
##########################################################

H <- 10  # number of events
N <- 1000  # number of trials per event
X <- cbind(1,rnorm(H))  # design matrix
qX <- ncol(X)
beta <- c(-3,0.15)  # coefficients
p <- expit(X%*%beta)  # probability
z <- rbinom(H,N,p)  # observed successes
sum(z)

##########################################################
### Fit model
##########################################################

source('~/Documents/git/GLM/binomial.logit.glm.mcmc.R')
priors <- list(sigma.beta=5)
tune <- list(beta=0.25)
start <- list(beta=coef(glm(cbind(z,N-z) ~ 0+X, family=binomial("logit"))))
out1 <- binomial.logit.glm.mcmc(z,N=N,X,priors,start,tune,n.mcmc=50000)

matplot(out1$beta, type="l", lty=1);abline(h=beta,col=1:qX,lty=3)