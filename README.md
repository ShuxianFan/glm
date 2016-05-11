# GLM

MCMC algorithms for implementing generalized linear models

#### Contents

##### Generalized linear model for binary data using a probit link
- **binary.probit.glm.pdf.R**: Model statement and full conditional distributions
- **binary.probit.glm.sim.R**: Simulate data and fit model
- *binary.probit.glm.mcmc.R*: MCMC algorithm for parameter estimation
	- **binary.probit.glm.mcmc.R**: MCMC algorithm for parameter estimation


- **binomial.logit.glm.mcmc.R**: MCMC algorithm for estimating parameters of GLM for binomially distributed data with logit link
- **binomial.logit.glm.sim.R**: Simulate data according to model specification of binomial.logit.glm.mcmc.R
- **binomial.logit.glm.pdf.R**: Model statement and full conditional distributions corresponding to binomial.logit.glm.mcmc.R

- **poisson.glm.mcmc.R**: MCMC algorithm for estimating parameters of Poisson GLM for count data
- **poisson.glm.sim.R**: Simulate data according to model specification of poisson.glm.mcmc.R
- **poisson.glm.pdf.R**: Model statement and full conditional distributions corresponding to poisson.glm.mcmc.R

- **nb.glm.mcmc.R**: MCMC algorithm for estimating parameters of negative binomial GLM for count data
- **nb.glm.sim.R**: Simulate data according to model specification of nb.glm.mcmc.R
- **nb.glm.pdf.R**: Model statement and full conditional distributions corresponding to nb.glm.mcmc.R
