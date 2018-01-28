# Feature filtration Script
# Input - features 
# Somnath Tagore

# Install packages
install.packages("coda")
install.packages("mcmc")

library(coda)
library(mcmc)

# Time calculation
normal_start_time <- proc.time()

# Input feature file
dat <- read.csv("bayes_test_data_1.csv", sep=",")

namebank <- dat[,3]
namebank

dat1 <- dat[,11]
dat1
length(dat1)
mean(dat1)
dat <- dat[,13]

dat

# Basic statistics calculation
meandata <- mean(dat1)
meandata
standarddeviationdata <- sd(dat1)
standarddeviationdata

# ABC-MCMC function
ABC_acceptance <- function(par){
  
  # prior to avoid negative standard deviation
  if (par[2] <= 0) return(F) 
  
  # stochastic model generates a sample for given par
  samples <- rnorm(mean(dat1), mean =par[1], sd = par[2])
  
  # comparison with the observed summary statistics
  diffmean <- abs(mean(samples) - meandata)
  print(diffmean)
  diffsd <- abs(sd(samples) - standarddeviationdata)
  print(diffsd)
  if((diffmean < 0.1) & (diffsd < 0.2)) return(T) else return(F)
}
mean(dat1)

# Plugging this in a standard metropolis hastings MCMC, 
# with the metropolis acceptance exchanged for the ABC acceptance
run_MCMC_ABC <- function(startvalue, iterations){
  
  chain = array(dim = c(iterations+1,2))
  chain[1,] = startvalue
  
  for (i in 1:iterations){
    
    # proposalfunction
    proposal = rnorm(2.26,mean = chain[i,], sd= c(0.7,0.7))
    
    if(ABC_acceptance(proposal)){
      chain[i+1,] = proposal
    }else{
      chain[i+1,] = chain[i,]
    }
  }
  summary(chain)
  return(mcmc(chain))
}

# Final posterior probability
posterior <- run_MCMC_ABC(c(0.0016,0.0012),300000)
plot(posterior)
