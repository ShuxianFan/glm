#
#
# Bayesian binomial generalized linear model
#
# Function name: binomial.glm.MCMC
#
# Author: Brian M. Brost
# Contact: bmbrost@gmail.com
#
# Last updated: 02 APR 2016
#
# Model statement:
#	z_i ~ Binom(N_i,p_i)
#	logit(p_i) = x_i%*%beta
#	beta ~ N(0,sigma.beta^2*I)
#
# Reference:
#
# Required R packages: 
#
# Inputs:
#
# z - vector of length n containing the number of successes observed during event i.
#	Elements in z correspond to rows in the design matrix X. Note that the value of 
#	z[1] corresponds to X[1,], z[2] corresponds to X[2,], etc.
# N - vector of length n containing the number of trials during event i
# X - design matrix of dimension n x qX containing covariates (plus
#	intercept) for which inference is desired
# priors - list of priors containing the following elements:
#	1. sigma.beta - standard deviation of normal prior on beta
# start - list of starting values containing the following elements:
#	1. beta - vector of starting values for coefficients
# tune - List of tuning parameters containing the following elements:
#	1. beta - tuning parameter for Metropolis-Hasting update on beta
# adapt - Switch to enable adapative tuning (TRUE/FALSE)
# n.mcmc - number of desired MCMC iterations
#
#

binomial.glm.mcmc <- function(z,N,X,priors,start,tune,adapt=TRUE,n.mcmc=1000){

	###
	###  Libraries and Subroutines
	###

	# library(mvtnorm)

	get.tune <- function(tune,keep,k,target=0.44){  # adaptive tuning
		# a <- min(0.01,1/sqrt(k))
		a <- min(0.025,1/sqrt(k))
		exp(ifelse(keep<target,log(tune)-a,log(tune)+a))
	}
	
	logit <- function(x){
		log(x/(1-x))
	}

	expit <- function(x){
		exp(x)/(1+exp(x))
	}

	###
	###  Setup variable, starting values, and priors
	###
# browser()
	qX <- ncol(X)  # number of betas
	
	beta <- start$beta
	p <- expit(X%*%beta)  # binomial probability

	mu.beta <- rep(0,qX)  # mean of normal prior on beta
	sigma.beta <- priors$sigma.beta
	
	# Sigma.beta <- solve(crossprod(X))
	# c <- sigma.beta/min(diag(Sigma.beta))
	# Sigma.beta <- c*Sigma.beta  # variance-covariance matrix of
		# # normal prior on beta

	###
	### Create receptacles for output
	###
  
	beta.save <- matrix(0,n.mcmc,qX)  # rsf coefficients
	
	keep <- list(beta=0)
	keep.tmp <- keep  # track MH accpetance rate for adaptive tuning
	Tb <- 50  # frequency of adaptive tuning

	###
	###  Begin MCMC loop
	###

	for(k in 1:n.mcmc){

		if(k%%1000==0) cat(k,"");flush.console()	

		if(adapt==TRUE & k%%Tb==0) {  # Adaptive tuning
			keep.tmp <- lapply(keep.tmp,function(x) x/Tb)
			tune$beta <- get.tune(tune$beta,keep.tmp$beta,k)
			keep.tmp <- lapply(keep.tmp,function(x) x*0)
	   	} 	

		###
		### Update beta
		### 
		
		beta.star <- rnorm(qX,beta,tune$beta)
		p.star <- expit(X%*%beta.star)  # binomial probability
  		mh.0 <- sum(dbinom(z,1000,p,log=TRUE))+
			sum(dnorm(beta,mu.beta,sigma.beta,log=TRUE))
		mh.star <- sum(dbinom(z,1000,p.star,log=TRUE))+
			sum(dnorm(beta.star,mu.beta,sigma.beta,log=TRUE))
		if(exp(mh.star-mh.0)>runif(1)){
			beta <- beta.star
			p <- p.star
			keep$beta <- keep$beta+1
			keep.tmp$beta <- keep.tmp$beta+1
		}

		###
		###  Save samples 
	    ###
# browser()
		beta.save[k,] <- beta
	}
	
	keep$beta <- keep$beta/n.mcmc
	cat(paste("\nbeta acceptance rate:",round(keep$beta,2))) 

	###
	### Write output
	###

	list(beta=beta.save,keep=keep,
		z=z,N=N,X=X,priors=priors,start=start,tune=tune,n.mcmc=n.mcmc)
}