###
### Simulate binary data and fit model using probit.reg.mcmc.R
###

rm(list=ls())

T <- 100  # number of observations

# Define covariates
time <- c(0,cumsum(rgamma(T-1,shape=1.1,scale=4.5)))  # time covariate
hr <- ((time)-24*floor(time/24))
day <- ceiling(time/24)
X <- cbind(1,day,hr)  # Design matrix
X[,-1] <- scale(X[,-1])
qX <- ncol(X)

beta <- c(-0.5,1.5,0.5)  # Coefficients on X

# Simulate data
p <- pnorm(X%*%beta)  # probability of 1
hist(p);summary(p)
y <- rbinom(T,1,p)  # haulout indicator variable: 1=hauled-out, 0=at-sea
table(y)

plot(time,y,ylim=range(c(y,X%*%beta)))
lines(time,X%*%beta,col=2)

# Fit model
source('~/Documents/git/GLM/probit.glm.mcmc.R', chdir = TRUE)
start <- list(beta=beta)
priors <- list(mu.beta=rep(0,qX),Sigma.beta=diag(qX)*100)
out1 <- probit.glm.mcmc(y,X,priors,start,1000)

matplot(out1$beta,type="l",lty=1);abline(h=beta,col=1:3,lty=2)
boxplot(pnorm(out1$u),col=8,outline=FALSE)
points(y,col=3)
out1$DIC
