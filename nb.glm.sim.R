
rm(list=ls())

##########################################################
### Simulate data for negative binomial GLM
##########################################################

n <- 1000  # number of obserations to simulate
X <- cbind(1,rnorm(n))  # design matrix
qX <- ncol(X)
beta <- c(-2,3)  # coefficients
mu <- exp(X%*%beta)  # intensity of Poisson process
hist(mu)
alpha <- 1 # dispersion parameter; var(z)=mu+mu^2/alpha, i.e.,  alpha=\infty yields Poisson
z <- rnbinom(n,size=alpha,mu=mu)  # observed counts
hist(z,breaks=20)


##########################################################
### Fit model
##########################################################

source('~/Documents/git/GLM/nb.glm.mcmc.R')
priors <- list(sigma.beta=5,a=1,b=0.1)
# hist(rgamma(1000,1,0.01),breaks=100)  # vague parameterization
tune <- list(beta=0.1,alpha=0.1)
start <- list(beta=coef(glm(z ~ 0+X, family=poisson())),alpha=alpha)
out1 <- nb.glm.mcmc(z,X,priors,start,tune,adapt=TRUE,n.mcmc=15000)

matplot(out1$beta, type="l", lty=1);abline(h=beta,col=1:qX,lty=3)
matplot(out1$alpha, type="l", lty=1);abline(h=alpha,col=1,lty=3)