###
### This is work in progress. This implementation is unintegrated,
### i.e., the random effects have not been marginalized out
###


#
#
# Bayesian negative binomial generalized linear model
#
# Function name: nb.glm.2.MCMC
#
# Author: Brian M. Brost
# Contact: bmbrost@gmail.com
#
# Last updated: 23 MAR 2016
#
# Model statement:
#	z_i ~ NB(lambda_i,alpha)
#	lambda_i = log(x_i%*%beta)*e_i
#	beta ~ N(0,sigma.beta^2*I)
# e_i ~ gamma(r,1/r)	
#r~gamma(a,b)
# 
#	note: E[alpha]=r^2 and Var[alpha]=r^3
#
# Reference:
#
# Required R packages: 
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
#	2. r - Shape and inverse rate parameter of gamma prior for alpha
# start - List of starting values containing the following elements:
#	1. beta - Vector of starting values for coefficients
# tune - list of tuning parameters containing the following elements:
#	1. beta - Tuning parameter for Metropolis-Hastings update on beta
# 2. r - Tuning parameter for Metropolis-Hastings update for random effects
# adapt - Switch to enable adapative tuning (TRUE/FALSE)
# n.mcmc - Number of desired MCMC iterations
#
#

nb.glm.2.mcmc <- function(z,X,priors,start,tune,adapt=TRUE,n.mcmc=1000){

	###
	###  Libraries and Subroutines
	###

	get.tune <- function(tune,keep,k,target=0.44){  # adaptive tuning
		# a <- min(0.01,1/sqrt(k))
		a <- min(0.025,1/sqrt(k))
		exp(ifelse(keep<target,log(tune)-a,log(tune)+a))
	}
	
	
	###
	###  Setup variable, starting values, and priors
	###

	qX <- ncol(X)  # number of betas
	n <- length(z)
	
	beta <- start$beta
	r <- start$r  # inverse dispersion parameter
	eps <- rgamma(n,shape=r,scale=1/r)
	lambda <- exp(X%*%beta)*eps  # mean of negative binomial
	
	mu.beta <- rep(0,qX)  # mean of normal prior on beta
	sigma.beta <- priors$sigma.beta
	a <- priors$a
	b <- priors$b
	
	# Sigma.beta <- solve(crossprod(X))
	# c <- sigma.beta/min(diag(Sigma.beta))
	# Sigma.beta <- c*Sigma.beta  # variance-covariance matrix of
		# # normal prior on beta

	###
	### Create receptacles for output
	###
  
	beta.save <- matrix(0,n.mcmc,qX)  # coefficients
	eps.save <- matrix(0,n.mcmc,n)  # random effects
	r.save <- numeric(n.mcmc)  # inverse dispersion parameter
	
	keep <- list(beta=0,r=0,eps=0)
	keep.tmp <- keep  # track MH accpetance rate for adaptive tuning
	Tb <- 50  # frequency of adaptive tuning

	###
	###  Begin MCMC loop
	###

	for(k in 1:n.mcmc){

		if(k%%1000==0) cat(k,"");flush.console()	

		if(adapt==TRUE & k%%Tb==0) {  # Adaptive tuning
			keep.tmp[1:2] <- lapply(keep.tmp[1:2],function(x) x/Tb)
			keep.tmp$eps <- keep.tmp$eps/(n*Tb)
			tune$beta <- get.tune(tune$beta,keep.tmp$beta,k)
			tune$r <- get.tune(tune$r,keep.tmp$r,k)
			tune$eps <- get.tune(tune$eps,keep.tmp$eps,k)
			keep.tmp <- lapply(keep.tmp,function(x) x*0)
	   	} 	

		###
		### Update beta
		### 
		
		beta.star <- rnorm(qX,beta,tune$beta)
		lambda.star <- exp(X%*%beta.star)*eps  # mean of negative binomial
		mh.0.beta <- sum(dpois(z,lambda,log=TRUE))+
		  sum(dnorm(beta,mu.beta,sigma.beta,log=TRUE))
		mh.star.beta <- sum(dpois(z,lambda.star,log=TRUE))+
		  sum(dnorm(beta.star,mu.beta,sigma.beta,log=TRUE))
		print(exp(mh.star.beta-mh.0.beta))
		if(exp(mh.star.beta-mh.0.beta)>runif(1)){
			beta <- beta.star
			lambda <- lambda.star
			keep$beta <- keep$beta+1
			keep.tmp$beta <- keep.tmp$beta+1
		}

		###
		### Update eps (random effect)
		### 
		
		eps.star <- rnorm(n,eps,tune$eps)
    idx1 <- which(eps.star>0)
    lambda.star <- exp(X[idx1,]%*%beta.star)*eps.star[idx1]  # mean of negative binomial
    mh.star.eps <- dpois(z[idx1],lambda.star,log=TRUE)+
		  dgamma(eps.star[idx1],shape=r,scale=1/r,log=TRUE)
		mh.0.eps <- dpois(z[idx1],lambda[idx1],log=TRUE)+
		  dgamma(eps[idx1],shape=r,scale=1/r,log=TRUE)
		idx2 <- which(exp(mh.star.eps-mh.0.eps)>runif(length(idx1)))
		eps[idx1[idx2]] <- eps.star[idx1[idx2]]
		lambda[idx1[idx2]] <- lambda.star[idx2]
		keep$eps <- keep$eps + length(idx2)
		
		###
		### Update r (inverse dispersion parameter)
		### 

		r.star <- rnorm(1,r,tune$r)
		if(r.star>0){
		  mh.0.r <- sum(dgamma(eps,shape=r,scale=1/r,log=TRUE))+
				dgamma(r,shape=a,rate=b,log=TRUE)
			mh.star.r <- sum(dgamma(eps,shape=r.star,scale=1/r.star,log=TRUE))+
				dgamma(r.star,shape=a,rate=b,log=TRUE)
			if(exp(mh.star.r-mh.0.r)>runif(1)){
				r <- r.star
				keep$r <- keep$r+1
				keep.tmp$r <- keep.tmp$r+1
			}		
		}

		###
		###  Save samples 
	  ###

		beta.save[k,] <- beta
		r.save[k] <- r
		eps.save[k,] <- eps
	}
	
	keep$beta <- keep$beta/n.mcmc
	keep$eps <- keep$eps/(n.mcmc*n)
	keep$r <- keep$r/n.mcmc
	cat(paste("\nbeta acceptance rate:",round(keep$beta,2))) 
	cat(paste("\neps acceptance rate:",round(keep$eps,2))) 
	cat(paste("\nr acceptance rate:",round(keep$r,2))) 
	
	###
	### Write output
	###

	list(beta=beta.save,r=r.save,eps=eps.save,keep=keep,
		z=z,X=X,priors=priors,start=start,tune=tune,n.mcmc=n.mcmc)
}