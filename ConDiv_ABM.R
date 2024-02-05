
###########################################################################################
## Deffner, D., Fedorova, N., Andrews, J. & McElreath, R.
## Bridging theory and data: A computational workflow for cultural evolution
## Agent-based simulation model (authored by D.Deffner; email: deffner@mpib-berlin.mpg.de)
###########################################################################################

#Load functions
library(scales)
library(RColorBrewer)
library(parallel)

#Define exponential function for age-dependent probabilities
#We assume that both survival and the probability to update variant decline as agents age
f_age <- function(increasing, rate, x) {
  if(increasing == TRUE){
    1-exp(-rate*(x-1))
  } else {
    exp(-rate*(x-1))
  }
}

#Define function to calculate 2D distance
Euclidian_distance <- function(x1,x2,y1,y2) sqrt( (x1-x2)^2 + (y1-y2)^2 )

#Define constant parameter values (i.e., parameters we do not perform a sweep for)
N = 3000         #Number of agents
N_groups = 30    #Number of groups
N_per_group = N/N_groups
N_burn_in = 1000 #Number of time steps to get equilibrium age structure
max_age <- 90     #Maximum Age
N_mod = 30       #Number of role models
m_NL = FALSE     #If TRUE, we take real migration rates from NL; if FALSE corresponding "effective" mig. rate
mu = 0.05        #Innovation rate (fraction of learning events that are innovations)
size_grid = 50   #Grid size for spatially-explicit model with local migration
r_dist  = 0      #Effect of distance on migration (if r_dist = 0, migrants choose village randomly)
r_learn = 0.03    #Decay rate for learning
r_mort  = 0.001   #Decay rate for survival

#Define parameter grid to loop over for sweep
seq<-expand.grid(N_steps=300, Nsim = 100, theta= seq(0.8, 5, 0.2), const_m = seq(0, 0.4, 0.02) )

#Define simulation function
sim.funct <- function(N_steps, Nsim, theta, const_m){
  
  #Overall output file
  Combined_list <- list()
  
  #Loop over independent simulations
  for (sim in 1:Nsim) {
    
    #Assign migration rates, depending on whether we use data from the NL
    if (m_NL == TRUE){
      m <- age_mig_NL
    } else {
      m <- rep(const_m, max_age)
    }
    
    #Create rectangular grid for villages, we try to make it spatially explicit
    grid <- matrix(0, nrow = size_grid, ncol = size_grid)
    
    # Now let's populate the grid with villages
    Locations <- sample(size_grid * size_grid, size = N_groups)
    grid[Locations] <- sample(N_groups)
    
    # Create (Euclidian) distance matrix between villages
    
    #Distance matrix
    dist <- matrix(0, nrow = N_groups, ncol = N_groups)
    for (x in 1:N_groups) {
      for (y in 1:N_groups) {
        dist[x,y] <- sqrt( (which(grid == x, arr.ind = TRUE)[1] - which(grid == y, arr.ind = TRUE)[1])^2 +
                             (which(grid == x, arr.ind = TRUE)[2] - which(grid == y, arr.ind = TRUE)[2])^2 )
      }
    }
    
    # Initialize population; we just keep track of personal id, group id and conformity exponent per individual
    Age <- sample(1:80, N, replace = TRUE)
    group <- rep(1:N_groups,each= N_per_group)
    
    #Initialize cultural traits 
    Traits <- rep(0, N)
    
    #Unique variant per group
    unique_per_group <- sample(1:N_groups,  replace = FALSE)
    for (g in 1:N_groups) Traits[group == g] <- unique_per_group[g]
    
    #Counter for cultural variants to make sure innovation actually produces new variant
    Counter <- max(unique_per_group)
    
    #Vector to store cultural Fst values
    Diversity <- rep(0,N_steps )
    
    #Loop over years (we first use N_burn_in time steps to reach demographic equilibrium)
    for (t in 1: (N_burn_in + N_steps ) ){
      print(t)
      # 1) Demographics
      # a) Birth-death process  
      for (g in 1:N_groups){
        idx_group <- which(group == g)
        
        #Create survivors
        alive <- rbinom(N_per_group , 1, f_age(FALSE, r_mort, Age[idx_group]) ) 
        
        #Limit life span to maximal age
        alive[which(Age[idx_group] == max_age)] <- 0
        
        #Select to fill empty slots
        babies  <- which(alive == 0)
        parents <- sample(which(alive == 1 & Age[idx_group] >= 18), length(babies), replace = TRUE)
        
        #Children copy traits of their parents
        Traits[idx_group[babies]] <- Traits[idx_group[parents] ]
        
        #Set children ages to 1 and increase rest by 1
        Age[idx_group[babies]] <- 1
        Age[idx_group[-babies]] <- Age[idx_group[-babies]] + 1
      }
      
      #If we're past the burn in, we include migration and cultural transmission
      if (t > N_burn_in){
        
        # b) Migration  
        #Create pool of migrants, according to (age-specific) migration rates
        migrate <- rbinom(N , 1, m[Age] ) 
        Migrants <- sample(which(migrate==1))
        
        #How many spots are available in each group
        Spots_per_group <- sapply(1:N_groups, function(i) length(which(group[Migrants] == i) ))
        
        #Vector to store new group id for each migrant
        new_group <- rep(0, length(Migrants))
        
        #If there are migrants, loop over them
        if (length(Migrants) > 0){
          for (i in Migrants)  {
            #Which groups are already full
            full <- which(Spots_per_group == 0)
            
            #If only one group is left, choose this one
            if (length(full) == (N_groups-1)){
              new <- which(Spots_per_group > 0)  
              
              #Otherwise, choose one group depending on distance  
            }else{
              probs <- exp(-r_dist * dist[group[i],])
              
              #Set prob for own group and for full groups to 0
              probs[full] <- 0
              probs[group[i]] <- 0
              probs <- probs/sum(probs)
              
              #Sample new group
              new <- sample((1:N_groups), 1, prob = probs )
            }
            
            #Assign new group for migrant and reduce number of free spots
            new_group[which(Migrants == i)] <- new
            Spots_per_group[new] <- Spots_per_group[new]-1
          }  
          
          #Assign new groups
          group[Migrants] <- new_group
        }
        
        # 2) Cultural Transmission
        
        #Create pool of learners
        Learners <- which( rbinom(N , 1, f_age(FALSE, r_learn, Age) ) == 1 )
        
        #Create vector for new variants
        Traits_new <- rep(0, length(Learners))
        
        #Loop over all learners
        for (i in Learners) {
          
          #Sample models
          Model_ids <- sample(which(group == group[i]), N_mod)
          
          #Vector with unique variants
          Variants <- unique(Traits[Model_ids]) 
          
          #Frequency of each variant
          Freq_Variants <- c()
          for (x in Variants) Freq_Variants[which(Variants == x)] <- length(which(Traits[Model_ids] == x))
          
          #Probability individuals choose each variant
          prob <- Freq_Variants^theta/ sum(Freq_Variants^theta)
          
          #Innovate new trait with probability mu
          if (runif(1)<mu){
            Traits_new[which(Learners == i)] <- Counter + 1
            
            #Update counter
            Counter <- Counter + 1
            
            #Socially learn with probability 1-mu    
          } else {
            if (length(Variants) == 1 ){
              Traits_new[which(Learners == i)] <- Variants
            } else {
              Traits_new[which(Learners == i)] <- sample(Variants, size = 1, prob = prob)
            }
          }
        }#i
        
        #Replace old by new cultural traits
        Traits[Learners] <- Traits_new
        
        #Quantify diversity (based on Mesoudi, 2018, Migration, acculturation, and the maintenance of between-group cultural variation)
        J <- unique(Traits)
        frequencies <- matrix(0, nrow = N_groups, ncol = length(J))
        for (g in 1:N_groups) frequencies[g,] <- sapply(J, function(j) length(which(Traits[which(group==g)] == j)) / N_per_group )
        
        total.var <- 1 - sum(colMeans(frequencies)^2)  # 1 - sum of squared means of each trait
        within.var <- mean(1 - rowSums(frequencies^2)) # mean of each group's (1 - the sum of squared freq of each trait)
        
        Diversity[t- N_burn_in] <- (total.var - within.var) / total.var
        
      }
      
    }#t
    
    Combined_list[[sim]]<- Diversity
    
  }#nsim
  
  return(Combined_list)  
  
}#end function  


#Pass to mclapply; it makes sense to select as many cores as there are parameter combinations in case you have access to a computer cluster ("mc.cores" argument)

result <- mclapply( 1:nrow(seq), function(i) sim.funct(seq$N_steps[i], seq$Nsim[i], seq$theta[i], seq$const_m[i]),mc.cores=77)


###
##
# Create plot (Fig. 3 in the manuscript)
##
###

#graphics.off()

pdf("Abstract_ABM.pdf", width = 11, height = 4)
#Crete color palette
col.pal <- brewer.pal(9, "Set1")

par( mar = c(4,4,0,0), oma = c(0,0,2.5,1))
layout(matrix(c(1,1,2,2,2,3,3), 1, 7, byrow = TRUE))

param_combi <- which(seq$theta == 1 & seq$const_m == 0)

t_plot <- 100

plot(result[[param_combi]][[1]][1:t_plot], type = "n", ylim = c(0,1), xlab = "", ylab = "")
for (i in 1:unique(seq$Nsim)) {
  lines(result[[param_combi]][[i]][1:t_plot], col = alpha(col.pal[1], alpha = 0.3))
}
text(70, 0.3, "no migration \n unbiased", cex = 1, col = col.pal[1])

param_combi <- 3
for (i in 1:unique(seq$Nsim)) {
  lines(result[[param_combi]][[i]][1:t_plot], col = alpha(col.pal[2], alpha = 0.3))
}
text(70, 0.9, "no migration \n weak conformity", cex = 1, col = col.pal[2])

#m=0.1
param_combi <- 113
for (i in 1:unique(seq$Nsim)) {
  lines(result[[param_combi]][[i]][1:t_plot], col = alpha(col.pal[3], alpha = 0.3))
}
text(70, 0.01, "migration \n weak conformity", cex = 1, col = col.pal[3])


param_combi <- 115
for (i in 1:unique(seq$Nsim)) {
  lines(result[[param_combi]][[i]][1:t_plot], col = alpha(col.pal[4], alpha = 0.3))
}
text(70, 0.5, "migration \n strong conformity", cex = 1, col = col.pal[4])

mtext(side = 1, line = 2.5, "Simulation Year")
mtext(side = 2, line = 2.5, "Cultural Fst")
mtext(side = 3, line = 0.5, "Simulation Dynamics" )
mtext('a', side=3, line=1, at=1)

# Remove first 100 timesteps and calculate mean for each parameter combination
MeanFst <- matrix(NA, nrow = nrow(seq), ncol = unique(seq$Nsim) )
for (i in 1:nrow(seq)){
  for (j in 1:unique(seq$Nsim)) {
    result[[i]][[j]] <- result[[i]][[j]]
    MeanFst[i,j] <- mean(result[[i]][[j]])
  }
}
OverallFst <- apply(MeanFst, 1, mean)

z <- matrix(NA, nrow = length(unique(seq$theta[which(seq$theta <= 3)])), ncol = length(unique(seq$const_m)) )
for (i in unique(seq$theta[which(seq$theta <= 3)])) {
  for (j in unique(seq$const_m)) {
    z[which(unique(seq$theta)==i),  which(unique(seq$const_m)==j)] <- mean(OverallFst[which(seq$theta==i & seq$const_m==j) ])
    
  }
}

image(1: 12, 1:21, z, col = alpha(rev(heat.colors(1000)), alpha = 1), zlim= c(0,max(z)),xaxt="n",yaxt = "n", xlab="", ylab="")
axis(side=1, at=seq(1,12,1), labels=unique(seq$theta[which(seq$theta <= 3)]) )
axis(side=2, at=seq(1,21,4), labels=unique(seq$const_m)[seq(1,21,4)]) 
rect(1.5, 0.5, 2.5, 21.5, density = NA, col = alpha("grey", alpha = 0), border = "black", lwd = 3)

mtext(side = 2, line = 2.5, "Migration rate" )
mtext(side = 3, line = 0.5, "Cultural Fst per parameter combination")

mtext(side = 1, line = 2.5, "Conformity exponent")
text(2, 10,"Unbiased transmission", cex = 1.5, srt = 90)
mtext('b', side=3, line=1, at=0)

#Causal effects
par( mar = c(4,4,0,0))
theta <- 1
low_mig <- which(seq$const_m == 0.1 & seq$theta == theta)
high_mig <- which(seq$const_m == 0.2 & seq$theta == theta)

low <- c()
for (i in low_mig) low <- c(low, unlist(result[[i]]) )

high <- c()
for (i in high_mig) high <- c(high, unlist(result[[i]]) )

contrast <- high-low
dens <- density(contrast)
x1 <- min(which(dens$x >= quantile(contrast, 0.05)))  
x2 <- max(which(dens$x <  quantile(contrast, 0.95)))
plot(dens, xlim = c(-0.4, 0), ylim = c(0,200), type="n", ann = FALSE, bty = "n", yaxt = "n")
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha(col.pal[1],alpha = 0.9), border = NA))

x1 <- min(which(dens$x >= quantile(contrast, 0)))  
x2 <- max(which(dens$x <  quantile(contrast, 1)))
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha(col.pal[1],alpha = 0.2), border = NA))

theta <- 1.6
low_mig <- which(seq$const_m == 0.1 & seq$theta == theta)
high_mig <- which(seq$const_m == 0.2 & seq$theta == theta)

low <- c()
for (i in low_mig) low <- c(low, unlist(result[[i]]) )

high <- c()
for (i in high_mig) high <- c(high, unlist(result[[i]]) )

contrast <- high-low
dens <- density(contrast)
x1 <- min(which(dens$x >= quantile(contrast, 0.05)))  
x2 <- max(which(dens$x <  quantile(contrast, 0.95)))
par(new = TRUE)
plot(dens, xlim = c(-0.4, 0), ylim = c(0,20), type="n", ann = FALSE, bty = "n", yaxt = "n")
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha(col.pal[2],alpha = 0.9), border = NA))

x1 <- min(which(dens$x >= quantile(contrast, 0)))  
x2 <- max(which(dens$x <  quantile(contrast, 1)))
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha(col.pal[2],alpha = 0.2), border = NA))

theta <- 2
low_mig <- which(seq$const_m == 0.1 & seq$theta == theta)
high_mig <- which(seq$const_m == 0.2 & seq$theta == theta)

low <- c()
for (i in low_mig) low <- c(low, unlist(result[[i]]) )

high <- c()
for (i in high_mig) high <- c(high, unlist(result[[i]]) )

contrast <- high-low
dens <- density(contrast)
x1 <- min(which(dens$x >= quantile(contrast, 0.05)))  
x2 <- max(which(dens$x <  quantile(contrast, 0.95)))
par(new = TRUE)
plot(dens, xlim = c(-0.4, 0), ylim = c(0,30), type="n", ann = FALSE, bty = "n", yaxt = "n")
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha(col.pal[3],alpha = 0.9), border = NA))

x1 <- min(which(dens$x >= quantile(contrast, 0)))  
x2 <- max(which(dens$x <  quantile(contrast, 1)))
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha(col.pal[3],alpha = 0.2), border = NA))

legend("top", title = expression("Conformity exp."), c("1","1.6","2"), col = c(col.pal[1],col.pal[2],col.pal[3]), cex = 1.2,bty = "n", lwd = 6, lty = 1)
mtext(side = 1, line = 2.5, "M -> CFst")
mtext(side = 3, line = 0.5, "Causal Effects" )
mtext('c', side=3, line=1, at=0)
abline(v = 0, lty = 2, col = "lightgrey")
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend(0.45, 0.9, c("1","", "0.9","", "0.8","", "0.7","","0.6","","0.5","","0.4","","0.3","","0.2","","0.1", "","0"), col= heat.colors(1000)[c(1,50,100,150,200,250, 300,350,400,450,500,550,600,650,700,750,800,850,900,950, 1000)], xpd = TRUE, inset = c(0, 0), bty = "n", pch=15,cex = 1.1, pt.cex = 3.2)

#dev.off()

