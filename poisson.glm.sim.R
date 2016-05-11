
rm(list=ls())

##########################################################
### Simulate data for Poisson GLM
##########################################################

n <- 1000  # number of obserations to simulate
X <- cbind(1,rnorm(n))  # design matrix
qX <- ncol(X)
beta <- c(-2,3)  # coefficients
lambda <- exp(X%*%beta)  # intensity of Poisson process
z <- rpois(n,lambda)  # observed counts

##########################################################
### Fit model
##########################################################

source('~/Documents/git/GLM/poisson.glm.mcmc.R')
priors <- list(sigma.beta=5)
tune <- list(beta=0.1)
start <- list(beta=coef(glm(z ~ 0+X, family=poisson())))
out1 <- poisson.glm.mcmc(z,X,priors,start,tune,n.mcmc=5000)

matplot(out1$beta, type="l", lty=1);abline(h=beta,col=1:qX,lty=3)


