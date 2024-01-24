
###########################################################################################
## Deffner, D., Fedorova, N., Andrews, J. & McElreath, R.
## Bridging theory and data: A computational workflow for cultural evolution
## Simulation of synthetic data and power analysis (authored by J.Andrews; email: jeffrey_andrews@eva.mpg.de)
###########################################################################################

#Set working directory to GitHub repo
setwd("~/GitHub/CulturalEvolutionWorkflow")

# Load the lattice package
library("lattice")
sigma <- seq(.01, 10, length.out = 10) # Create a vector that represents the amount of unexplained variance in the diversity metric
bm <- 10^-(0:9) # Create a vector that varies the magnitude of the effect of migration (M) on diversity
sweep <- expand.grid(bm, sigma) # Create a Matrix for the parameter sweep
n <- 30 # Define the number of countries. 

# Make a function to simulate the fake data
sim_dat <-  function(bm, sigma){
            C <- rnorm(n, 0,  1) # Generate a standardized distribution for the amount of conformity in each group - this is data in the emprical case. 
            M <- rnorm(n, 0,  1) # Generate a standardized distribution for the amount of conformity in each group - this is data in the emprical case. 
            D <- rnorm(n, M*bm + C*1 , sigma) # Simulate the values of Diversity 
            return(data.frame(C, M, D))
}


##########################
###### POWER ANALYSIS ####
m <- matrix(NA, nrow = 10, ncol = 10) # intialize a matrix for storing the data from the power analysis           
for(i in 1:nrow(sweep)){
  print(i)
  data <- list() # initalize a list for storing data
  for(j in 1:200) data[[j]] <- sim_dat(sweep[i, 1], sweep[i, 2]) # simulate the data
  est <- c() # intalize a vector for storing the data
  for(j in 1:length(data))  est <- c(est, log(abs(coef(summary(lm(D ~ M + C, data[[i]])))["M", "Estimate"]/sweep[i, 1]))) # Calculate the log error in the estimates. 
  m[i] <- mean(est) # store the error 
}

# Plot the error
levelplot(m, xlab = "Sigma", ylab ="10^-(bm)", col.regions=rev(heat.colors(100)), main = "log(bm_real-bm_est)")

