probit.reg.mcmc <- function(y,X,priors,start,n.mcmc){

	###
	### Brian M. Brost (10 AUG 2015)
	### Regression for binary data using probit link
	###

	###
	### Model statement:
	### y=0,u_t<=0
	### y=1,u_t>0
	### u_t~N(x*beta,1)
	### beta~N(0,sigma.beta^2*I)
	###	

	###
	###  Libraries and Subroutines
	###
	
	truncnormsamp <- function(mu,sig2,low,high,nsamp){
	  flow=pnorm(low,mu,sqrt(sig2)) 
	  fhigh=pnorm(high,mu,sqrt(sig2)) 
	  u=runif(nsamp) 
	  tmp=flow+u*(fhigh-flow)
	  x=qnorm(tmp,mu,sqrt(sig2))
	  x
	}
	
	###
	###  Preliminary Variables
	###

	# X=as.matrix(X)
	# y=as.vector(y)
	n <- length(y)
	qX <- ncol(X)
	y1 <- (y==1)
	y0 <- (y==0)
	y1.sum <- sum(y1)
	y0.sum <- sum(y0)
	u <- numeric(n)

	###
	### Starting values and priors
 	###

	beta <- start$beta
	# beta <- matrix(glm(y~data.frame(X[,-1]),
		# family=binomial(link="probit"))$coefficients,qX,1)		
	mu.beta <- priors$mu.beta
	Sigma.beta <- priors$Sigma.beta
	Sigma.beta.inv <- solve(Sigma.beta)

	###
	### Create receptacles for output
	###
  
	beta.save <- matrix(0,n.mcmc,qX)
	u.save <- matrix(0,n.mcmc,n)
	D.bar.save <- numeric(n.mcmc)  # D.bar for DIC calculation
	
	###
	###  Gibbs Loop
	###
	
	for(k in 1:n.mcmc){
		if(k%%1000==0) cat(k," ");flush.console()
	
		###
		### Sample u (auxilliary variable for probit regression)
	  	###
	 
		linpred <- X%*%beta
	  	u[y1] <- truncnormsamp(linpred[y1],1,0,Inf,y1.sum)
	  	u[y0] <- truncnormsamp(linpred[y0],1,-Inf,0,y0.sum)
	
	  	###
	  	### Sample beta
	  	###
# browser()		
	  	A.inv <- solve(t(X)%*%X+Sigma.beta.inv)
	  	b <- t(X)%*%u  # +mu.beta%*%Sigma.beta.inv
	  	beta <- A.inv%*%b+t(chol(A.inv))%*%matrix(rnorm(qX),qX,1)

		# tmpvar <- solve(t(X)%*%X + solve(Sigma.beta))
	  	# tmpmean <- tmpvar%*%(t(X)%*%u + solve(Sigma.beta)%*%mu.beta)
	  	# beta <- tmpmean+t(chol(tmpvar))%*%matrix(rnorm(qX),qX,1)
	
	  	###
	  	### Save Samples 
	  	###
	
	  	beta.save[k,] <- beta
	  	u.save[k,] <- u
	  	D.bar.save[k] <- -2*(sum(dbinom(y,1,pnorm(u),log=TRUE)))
	}
		
	###
	###  Calculate DIC 
	###
	
	if(qX==1)  postbetamn <- mean(beta.save)
	if(qX>1)  postbetamn <- apply(beta.save,2,mean)
	postumn <- apply(u.save,2,mean)
	D.hat=-2*(sum(dbinom(y,1,pnorm(postumn),log=TRUE)))
	D.bar <- mean(D.bar.save)
	pD <- D.bar-D.hat
	DIC <- D.hat+2*pD
	 
	###
	###  Write output 
	###
	
	list(beta=beta.save,u=u.save,DIC=DIC,n.mcmc=n.mcmc)
}
