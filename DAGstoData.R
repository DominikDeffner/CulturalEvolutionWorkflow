
###########################################################################################
## Deffner, D., Fedorova, N., Andrews, J. & McElreath, R.
## Bridging theory and data: A computational workflow for cultural evolution
## Simulation of synthetic data and power analysis (authored by J.Andrews; email: jeffrey_andrews@eva.mpg.de)
###########################################################################################

#Set working directory to GitHub repo
setwd("~/GitHub/CulturalEvolutionWorkflow")

# Load the lattice package
library("lattice")
sigma <- seq(.01, 3, length.out = 10) # Create a vector that represents the amount of unexplained variance in the diversity metric
bmd <- 10^-(0:9) # Create a vector that varies the magnitude of the effect of migration (M) on diversity
bcd <- 1 #Effect of conformity on diversity
bam <- 1 #Effect of age on migration rates
bac <- 1 #Effect of age on conformity
  
sweep <- expand.grid(bmd, sigma) # Create a Matrix for the parameter sweep
n <- 30 # Define the number of countries. 

# Make a function to simulate the fake data
sim_dat <-  function(bmd, sigma){
            A <- rnorm(n, 0,  1) # Generate a standardized distribution for age - this is data in the empirical case. 
            C <- rnorm(n, A*bac,  1) # Generate a standardized distribution for the amount of conformity in each group - this is data in the empirical case. 
            M <- rnorm(n, A*bam,  1) # Generate a standardized distribution for the amount of migration in each group - this is data in the empirical case. 
            D <- rnorm(n, M*bmd + C*bcd , sigma) # Simulate the values of Diversity 
            return(data.frame(A, C, M, D))
}


##########################
###### POWER ANALYSIS ####
m <- matrix(NA, nrow = 10, ncol = 10) # intialize a matrix for storing the data from the power analysis   

for(i in 1:nrow(sweep)){
  print(i)
  data <- list() # initalize a list for storing data
  for(j in 1:200) data[[j]] <- sim_dat(sweep[i, 1], sweep[i, 2]) # simulate the data
  error <- c() # intalize a vector for storing the data
  est <- c() # intalize a vector for storing the data
  
  for(j in 1:length(data)){
    estimate <- coef(summary(lm(D ~ M + C, data[[i]])))["M", "Estimate"]
    est <- c(est, estimate)
    error <- c(error, log(abs(estimate/sweep[i, 1] ) )) # Calculate the log error in the estimates. 
  }  
  m[i] <- mean(error) # store the error 
}


# Plot the error
levelplot(m, xlab = "Unexplained variance (Sigma)", ylab ="Effect of migration [10^-(bmd)]", col.regions=rev(heat.colors(100)), main = "Error [log(bmd_real-bmd_est)]")

