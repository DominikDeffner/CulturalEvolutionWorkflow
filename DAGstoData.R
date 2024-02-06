
###########################################################################################
## Deffner, D., Fedorova, N., Andrews, J. & McElreath, R.
## Bridging theory and data: A computational workflow for cultural evolution
## Simulation of synthetic data and power analysis
## authored by J.Andrews (email: jeffrey_andrews@eva.mpg.de) and Dominik Deffner (deffner@mpib-berlin.mpg.de)
###########################################################################################

#Set working directory to GitHub repo
setwd("~/GitHub/CulturalEvolutionWorkflow")

# Load the lattice package
library("lattice")
library(rethinking)
library(scales)
library(RColorBrewer)

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
##### Here, we use maximum likelihood methods with lm() to speed up process

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



#########################
###### CAUSAL EFFECT ####
##### Here, we use MCMC with stan through the rethinking package

dat <- sim_dat(1, 1)

m <- ulam(
  alist(
    D ~ dnorm( mu , sigma ) ,
    mu <- a + bmd*M + bcd*C,
    c(a,bmd,bcd) ~ dnorm( 0 , 1 ) ,
    sigma ~ dexp( 1 )
  ) , data=dat , chains=4, iter = 3000 )
s <- extract_post_ulam(m)

#Compute causal effect of migration on diversity
N = length(s$sigma)
Causal_effect <- matrix(NA, 3,N)
for (i in 1:3) {
  C <- c(-1, 0, 1)[i]
  Causal_effect[i,] <- rnorm(N, s$a+ s$bcd*C + s$bmd * 1, s$sigma) - rnorm(N, s$a+ s$bcd*C + s$bmd * (-1), s$sigma)
}

par(mar = c(2,0,0,0), oma = rep(1,4))
col.pal <- brewer.pal(9, "Set1")
dens <- density(Causal_effect[1,])
x1 <- min(which(dens$x >= quantile(Causal_effect[1,], 0.05)))  
x2 <- max(which(dens$x <  quantile(Causal_effect[1,], 0.95)))
plot(dens, xlim = c(-4, 8), ylim = c(0,0.4), type="n", ann = FALSE, bty = "n", yaxt = "n")
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha(col.pal[1],alpha = 0.9), border = NA))

x1 <- min(which(dens$x >= quantile(Causal_effect[1,], 0)))  
x2 <- max(which(dens$x <  quantile(Causal_effect[1,], 1)))
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha(col.pal[1],alpha = 0.2), border = NA))

dens <- density(Causal_effect[2,])
x1 <- min(which(dens$x >= quantile(Causal_effect[2,], 0.05)))  
x2 <- max(which(dens$x <  quantile(Causal_effect[2,], 0.95)))
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha(col.pal[2],alpha = 0.9), border = NA))

x1 <- min(which(dens$x >= quantile(Causal_effect[2,], 0)))  
x2 <- max(which(dens$x <  quantile(Causal_effect[2,], 1)))
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha(col.pal[2],alpha = 0.2), border = NA))

dens <- density(Causal_effect[3,])
x1 <- min(which(dens$x >= quantile(Causal_effect[3,], 0.05)))  
x2 <- max(which(dens$x <  quantile(Causal_effect[3,], 0.95)))
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha(col.pal[3],alpha = 0.9), border = NA))

x1 <- min(which(dens$x >= quantile(Causal_effect[3,], 0)))  
x2 <- max(which(dens$x <  quantile(Causal_effect[3,], 1)))
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha(col.pal[3],alpha = 0.2), border = NA))
abline(v = 0, lty = 2, col = "lightgrey")
legend("topleft", title = expression("Level of Conformity"), c("-1","0","1"), col = c(col.pal[1],col.pal[2],col.pal[3]), cex = 1.2, bty = "n", lwd = 6, lty = 1)
mtext(side = 1, line = 2, "M -> D", cex = 1.2)
mtext(side = 3, line = 1, "Causal Effects", cex = 1.2)
