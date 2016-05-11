#
#
# Bayesian Poisson generalized linear model for count data
#
# Function name: poisson.glm.MCMC
#
# Author: Brian M. Brost
# Contact: bmbrost@gmail.com
#
# Last updated: 18 MAR 2016
#
# Model statement:
#	z_i ~ Pois(lambda_i)
#	log(lambda_i) = x_i%*%beta
#	beta ~ N(0,sigma.beta^2*I)
#
# Reference:
#
# Required R packages: mvtnorm (if using multivariate normal prior on beta)
#
# Inputs:
#
# z - Vector of length n containing the count of observations corresponding to 
# 	each row in the design matrix X. Note that the value of z[1] corresponds to
#	X[1,], z[2] corresponds to X[2,], etc.
# X - Design matrix of dimension n x qX containing covariates (plus
#	intercept) for which inference is desired
# priors - List of priors containing the following elements:
#	1. sigma.beta - Standard deviation of normal prior on beta
# start - List of starting values containing the following elements:
#	1. beta - Vector of starting values for coefficients
# tune - List of tuning parameters containing the following elements:
#	1. beta - Tuning parameter for Metropolis-Hasting update on beta
# adapt - Switch to enable adapative tuning (TRUE/FALSE)
# n.mcmc - Number of desired MCMC iterations
#
#

poisson.glm.mcmc <- function(z,X,priors,start,tune,adapt=TRUE,n.mcmc=1000){

	###
	###  Libraries and Subroutines
	###

	# library(mvtnorm)

	get.tune <- function(tune,keep,k,target=0.44){  # adaptive tuning
		# a <- min(0.01,1/sqrt(k))
		a <- min(0.025,1/sqrt(k))
		exp(ifelse(keep<target,log(tune)-a,log(tune)+a))
	}
	
	###
	###  Setup variable, starting values, and priors
	###

	qX <- ncol(X)  # number of betas
	
	beta <- start$beta
	lambda <- exp(X%*%beta)  # intensity of Poisson process

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
		lambda.star <- exp(X%*%beta.star)  # intensity of Poisson process
  		mh.0 <- sum(dpois(z,lambda,log=TRUE))+
			sum(dnorm(beta,mu.beta,sigma.beta,log=TRUE))
		mh.star <- sum(dpois(z,lambda.star,log=TRUE))+
			sum(dnorm(beta.star,mu.beta,sigma.beta,log=TRUE))
		if(exp(mh.star-mh.0)>runif(1)){
			beta <- beta.star
			lambda <- lambda.star
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
		z=z,X=X,priors=priors,start=start,tune=tune,n.mcmc=n.mcmc)
}