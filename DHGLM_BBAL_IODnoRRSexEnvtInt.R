
model {
  #priors
  ## fixed effects
  for (j in 1:6) {
    beta[j] ~ dnorm(0, 0.001)     #for the mean part. beta[1] being the intercept now
  }
  
  for (k in 1:6) {
    gamma[k] ~ dnorm(0, 0.001)     #for the dispersion part
  }

### new estimates for Sexe interaction with the environment. Now we have an estimate environmental slope for each sexe
  for (j in 1) {			
    betaSi[j] ~ dnorm(0, 0.001)		# interaction environment * sexe
    gammaSi[j] ~ dnorm(0, 0.001)
  }

### new estimates for reproductive status effects interaction with the environment. Now we have an estimate of environmental slope for each reproductive status 
  for (j in 1:2) {
    betaRi[j] ~ dnorm(0, 0.001) 
    gammaRi[j] ~ dnorm(0, 0.001) 
  }

   ##random effects
  #you need to define a prior for random effect including the distribution of the variance and a sampling of the data for each group/individual
  ### Cycle (years)
 # updated priors for pair effect included now both in mean and variance part
  tau.yr[1:2, 1:2] ~ dwish(mat.yr[, ], 3) #wishart prior on the precision matrix
  sig2.yr[1:2, 1:2] <- inverse(tau.yr[, ]) #conversion to a variance matrix
  for (l in 1:n.Cycles) {
     uyr[l, 1:2] ~ dmnorm(zeroyr[], tau.yr[, ])
  }
    
  ###PairBand
 # updated priors for pair effect included now both in mean and variance part
  tau.pair[1:2, 1:2] ~ dwish(mat.pair[, ], 3) #wishart prior on the precision matrix
  sig2.pair[1:2, 1:2] <- inverse(tau.pair[, ]) #conversion to a variance matrix
  for (l in 1:n.Pairs) {
     upair[l, 1:2] ~ dmnorm(zerop[], tau.pair[, ])
  }
  
  
   
  ###indivual identity
  tau.id[1:2, 1:2] ~ dwish(mat.id[, ], 3) #wishart prior on the precision matrix
  sig2.id[1:2, 1:2] <- inverse(tau.id[, ]) #conversion to a variance matrix
  for (o in 1:n.ids) {
  uid[o, 1:2] ~ dmnorm(zero[], tau.id[, ]) #values for each individual. This is now a matrix (2 values per individuals)
  }
  
  #likelihood(i.e. model)
  for (i in 1:n) {
    #define the model for a gaussian variable with a mean and variance to estimate
    y[i] ~ dnorm(y.hat[i], tau.y[i])
    
    #model for the mean y.hat
    y.hat[i] <-
      beta[1] + beta[2] * Sexe[i] + beta[3] * Age[i] + beta[4] * Brooding[i] + beta[5] * Preponte[i]
      + (beta[6] + betaSi[1] * Sexe[i] + betaRi[1] * Brooding[i] + betaRi[2] * Preponte[i] ) * IODAnnual[i] 
        + uid[id[i], 1] + upair[PairBand[i], 1] + uyr[Cycle[i], 1]
    
    #model for the variance
    tau.y[i] <- 1 / sigma2[i] #expressing the model as variance instead of a precision

    log(sigma2[i])  <- gamma[1] + gamma[2] * Sexe[i] + gamma[3] * Age[i] + gamma[4] * Brooding[i] + gamma[5] * Preponte[i]
	+ (gamma[6] + gammaSi[1] * Sexe[i] + gammaRi[1] * Brooding[i] + gammaRi[2] * Preponte[i]) * IODAnnual[i]
	+ uid[id[i], 2] + upair[PairBand[i], 2] + uyr[Cycle[i], 2]
  }
   
}
 
