######################################
#
# COVID Monitoring
# Case Reporting w/ Delay
#
# County Model
#
# 07/28/21
#
######################################

library(nimble)
library(coda)


set.seed(576476)

#####################################
#Functions
#####################################
# Beta-Binomial distribution functions.

dbetabin=nimbleFunction(run=function(x=double(0),mu=double(0),phi=double(0),size=double(0),log=integer(0)){
  returnType(double(0))
  if(x>=0&x<=size){
    return(lgamma(size+1)+lgamma(x+mu*phi)+lgamma(size-x+(1-mu)*phi)+lgamma(phi)-
             lgamma(size+phi)-lgamma(mu*phi)-lgamma((1-mu)*phi)-lgamma(size-x+1)-lgamma(x+1))
  }else{
    return(-Inf)
  }
})

rbetabin=nimbleFunction(run=function(n=integer(0),mu=double(0),phi=double(0),size=double(0)){
  pi=rbeta(1,mu*phi,(1-mu)*phi)
  returnType(double(0))
  return(rbinom(1,size,pi))
})



#Replace MONTH with desired month
load('data_MONTH.Rda')

#Rda file contains all files needed to execute MCMC in nimble

#Model Code contains the model specification

#Constants contains:
#N - 91 (90 day window + 1 for next day forecast
#C - 89 (89 days with partial data, data for day 90 ignored and forecasted) 
#D - 30 (maximum reporting delay) 
#L - 88 (number of counties) 
#num - number of neighbors for each county
#adj - vector defining adjacency

#Model data contains:
#z - reporting matrix (dimensions - onset date, reporting delay, county) 
#y - case time series (dimensions - onset date, county) 
#off - log(population) for each county
#X - design matrix for day of the week

#Inits contains initial values for the MCMC


#Build the model.
model=nimbleModel(model_code,constants,model_data,inits)

#Compile the model
compiled_model=compileNimble(model,resetFunctions = TRUE)

#Set monitors
mcmc_conf=configureMCMC(model,monitors=c('lambda','alpha','delta','y','tau.dc','d0','d.c'),useConjugacy = TRUE)

#Build MCMC
mcmc<-buildMCMC(mcmc_conf)

#Compile MCMC
compiled_mcmc<-compileNimble(mcmc, project = model)

#Run the model 
samples=runMCMC(compiled_mcmc,inits=inits,
                nchains = 1, nburnin=15000,niter = 30000,samplesAsCodaMCMC = TRUE,thin=10,
                summary = FALSE, WAIC = FALSE,progressBar=TRUE) 

