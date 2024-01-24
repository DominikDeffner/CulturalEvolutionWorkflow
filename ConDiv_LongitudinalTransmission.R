
###########################################################################################
## Deffner, D., Fedorova, N., Andrews, J. & McElreath, R.
## Bridging theory and data: A computational workflow for cultural evolution
## Longitudinal transmission analysis (authored by D.Deffner; email: deffner@mpib-berlin.mpg.de)
###########################################################################################

#This script simulates longitudinal data on cultural traits and social networks for hypothetical participants
#It then fits a time-series transmission model in stan and plots the results (Fig.4 in the manuscript)

#Set working directory to GitHub repo
setwd("~/GitHub/CulturalEvolutionWorkflow")

#Load functions
library(scales)
library(RColorBrewer)

#Load real age trajectory for migration (based on Fedorova et al., 2022, The complex life course of mobility: Quantitative description of 300,000 residential moves in 1850-1950 Netherlands) ) 
age_mig_NL <- readRDS("beta_df.RDS")$beta_mod_f
max_age <- length(age_mig_NL)

#Define exponential function for age-dependent probabilities
#We assume that both survival and the probability to update your variant decline as agents age
f_age <- function(increasing, rate, x) {
  if(increasing == TRUE){
    1-exp(-rate*(x-1))
  } else {
    exp(-rate*(x-1))
  }
}

#Function for getting fixation Index (Fst) from trait vectors (based on Mesoudi, 2018, Migration, acculturation, and the maintenance of between-group cultural variation)
getFst <- function(trait_vec) {
  J <- unique(Traits)
  frequencies <- matrix(0, nrow = N_groups, ncol = length(J))
  for (g in 1:N_groups) frequencies[g,] <- sapply(J, function(j) length(which(trait_vec[which(group==g)] == j)) / N_per_group )
  
  total.var <- 1 - sum(colMeans(frequencies)^2)  # 1 - sum of squared means of each trait
  within.var <- mean(1 - rowSums(frequencies^2)) # mean of each group's (1 - the sum of squared freq of each trait)
  (total.var - within.var) / total.var  # return Fst
}



N = 3000         #Number of agents
N_groups = 30    #Number of groups
N_per_group = N/N_groups
N_steps = 30    #Number of time steps
N_burn_in = 1000 #Number of time steps to get equilibrium age structure
N_skip = 100     #Number of time steps before we record choices for analysis
N_mod = 10       #Number of role models
m_NL = TRUE      #If TRUE, we take real migration rates from NL; if FALSE corresponding "effective" mig. rate
m_const = 0.01
mu = 0.1        #Innovation rate (fraction of learning events that are innovations)
f = 3            #Conformity parameter
size_grid = 50   #Grid size for spatially-explicit model with local migration
r_dist  = 0      #Effect of distance on migration (if r_dist = 0, migrants choose village randomly)

#Define developmental rates for learning and mortality
r_learn = 0
r_mort  = 0.001


if (m_NL == TRUE){
  m <- age_mig_NL
} else {
  # m <- rep( sum(age_mig_NL*Age_hat), max_age )
  m <- rep(m_const, max_age)
}


#First create rectangular grid for villages, we try to make it spatially explicit
grid <- matrix(0, nrow = size_grid, ncol = size_grid)

# Now let's populate the grid with villages
Locations <- sample(size_grid * size_grid, size = N_groups)
grid[Locations] <- sample(N_groups)

# Create (Euclidian) distance matrix between villages
Euclidian_distance <- function(x1,x2,y1,y2) sqrt( (x1-x2)^2 + (y1-y2)^2 )

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
id <- sample(1:N)

#Initialize cultural traits 
Traits <- rep(0, N)

#Unique
unique_per_group <- sample(1:N_groups,  replace = FALSE)
for (g in 1:N_groups) Traits[group == g] <- unique_per_group[g]

Counter <- max(unique_per_group)


#Initialize data objects
#How many individuals to record data from

N_record = 200

#Traits
dat_trait  <- matrix(NA, nrow = N_record, ncol = N_steps)

#Interaction partners
dat_models <- array(NA, dim = c(N_record, N_steps, N_mod))

#Ages
dat_age  <- matrix(NA, nrow = N_record, ncol = N_steps)

#Group
dat_group  <- matrix(NA, nrow = N_record, ncol = N_steps)

#Initialize Fst vector
Diversity <- rep(0,N_steps )


#Loop over years
for (t in 1: (N_burn_in + N_skip + N_steps ) ){
  print(t)
  
  # 1) Demographics
  # a) Birth-death process  
  for (g in 1:N_groups){
    idx_group <- which(group == g)
    
    #Create survivors
    alive <- rbinom(N_per_group , 1, f_age(FALSE, r_mort, Age[idx_group]) ) 
    
    #Limit life span to maximal age from data
    alive[which(Age[idx_group] == max_age)] <- 0
    
    #Select to fill empty slots
    babies  <- which(alive == 0)
    parents <- sample(which(alive == 1 & Age[idx_group] >= 18), length(babies), replace = TRUE)
    
    #Children copy traits of their parents
    Traits[idx_group[babies]] <- Traits[idx_group[parents] ]
    
    #Update ids for newborns
    if (t > (N_burn_in+N_skip)){
      id[idx_group[babies]] <- (max(id)+1):(max(id)+length(babies))
    } 
    
    
    #Set childrens' ages to 1 and increase rest by 1
    Age[idx_group[babies]] <- 1
    Age[idx_group[-babies]] <- Age[idx_group[-babies]] + 1
  }
  
  if (t > N_burn_in){
    
    # b) Migration  
    
    #Create pool of migrants, accoring to age-specific migration rates
    migrate <- rbinom(N , 1, m[Age] ) 
    Migrants <- sample(which(migrate==1))
    
    #How many spots are available in each group
    Spots_per_group <- sapply(1:N_groups, function(i) length(which(group[Migrants] == i) ))
    
    #Vector to store new group id for each migrant
    new_group <- rep(0, length(Migrants))
    
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
    
    #Loop over all individuals
    
    for (i in 1:N) {
      
      #Sample models
      Model_ids <- sample(which(group == group[i]), N_mod)
      
      if (t > (N_burn_in+N_skip)){
        if (id[i] %in% 1:N_record){
          dat_models[id[i], t - (N_burn_in+N_skip), ] <- Traits[Model_ids]
        } 
      }
      
      if (i %in% Learners){
        
        #Vector with unique variants including own variant
        Variants <- unique(Traits[c(Model_ids, i)]) 
        
        #Frequency of each variant
        Freq_Variants <- c()
        for (x in Variants) Freq_Variants[which(Variants == x)] <- length(which(Traits[c(Model_ids,i)] == x))
        
        #Probability individuals choose each variant
        prob <- Freq_Variants^f / sum(Freq_Variants^f)
        
        #Innovate new trait with probability mu
        if (runif(1)<mu){
          Traits_new[which(Learners == i)] <- Counter+1
          
          #Update Counter
          Counter <- Counter+1
          
          #Socially learn with probability 1-mu    
        } else {
          
          if (length(Variants) == 1 ){
            Traits_new[which(Learners == i)] <- Variants
          } else {
            Traits_new[which(Learners == i)] <- sample(Variants, size = 1, prob = prob)
          }
        }
        
      }#if i
      
    }#i
    
    #Replace old by new cultural traits
    Traits[Learners] <- Traits_new
    
    #Record traits
    
    if (t > (N_burn_in+N_skip)){
      #Quantify diversity (based on Mesoudi, 2018, Migration, acculturation, and the maintenance of between-group cultural variation)
      J <- unique(Traits)
      frequencies <- matrix(0, nrow = N_groups, ncol = length(J))
      for (g in 1:N_groups) frequencies[g,] <- sapply(J, function(j) length(which(Traits[which(group==g)] == j)) / N_per_group )
      
      total.var <- 1 - sum(colMeans(frequencies)^2)  # 1 - sum of squared means of each trait
      within.var <- mean(1 - rowSums(frequencies^2)) # mean of each group's (1 - the sum of squared freq of each trait)
      
      Diversity[t - (N_burn_in+N_skip)] <- (total.var - within.var) / total.var
      
      for (i in 1:N_record) {
        if (i %in% id){
          dat_trait[i, t - (N_burn_in+N_skip)] <- Traits[which(id == i)]
          dat_age[i, t - (N_burn_in+N_skip)]   <- Age[which(id == i)]
          dat_group[i, t - (N_burn_in+N_skip)] <- group[which(id == i)]
          
        } 
      }
    }
    
    
  }
  
  
  
}#t



#Restructure data for stan

choices  <- rep(-10, N_record*N_steps)
innovate <- rep(-10, N_record*N_steps)
id <- rep(1:N_record, each = N_steps)
group <- rep(-10, N_record*N_steps)
age <- rep(-10, N_record*N_steps)
migrate <- rep(-10, N_record*N_steps)
alternatives <- matrix(-10, nrow = N_record*N_steps, ncol = N_mod+1)
N_alt <- rep(-10, N_record*N_steps)
keep <- rep(-10, N_record*N_steps)

#Loop over individuals
for (i in 1:N_record) {
  for (t in 2:N_steps) {
    
    keep[which(id==i)][t] <- ifelse(is.na(dat_trait[i, t]), 0, 1)
    group[which(id==i)][t] <- dat_group[i, t]
    migrate[which(id==i)][t] <-ifelse(group[which(id==i)][t-1] == group[which(id==i)][t], 0, 1)
    
    age[which(id==i)][t] <- dat_age[i, t]
    
    choice <- ifelse(is.na(dat_trait[i, t]), -10, dat_trait[i, t])
    all_traits <- c(dat_models[i,t, ], dat_trait[i,t-1])
    
    available <- unique(all_traits)
    n <- sapply(available, function(x) length(which(all_traits == x)))
    
    alternatives[which(id==i)[t], 1:length(n)] <- n
    N_alt[which(id==i)][t] <- length(n)
    
    if (choice %in% available){
      innovate[which(id==i)][t] <- 0
      choices[which(id==i)][t] <- which(available == choice)
    }
    else if (choice == -10){
      innovate[which(id==i)][t] <- -10
    } else{
      innovate[which(id==i)][t] <- 1
      
    }
    
  }
}




#Stan data

stan.data <- list(choices = choices[keep == 1],
                  innovate  = innovate[keep == 1],
                  migrate  = migrate[keep == 1],
                  id  = id[keep == 1],
                  age  = age[keep == 1],
                  group  = group[keep == 1],
                  frequencies = alternatives[keep == 1,],
                  N_alt  = N_alt[keep == 1])

stan.data$N    <- length(stan.data$choices)
stan.data$N_partners  <- N_mod
stan.data$Max_age  <- max(stan.data$age)

ids <- unique(stan.data$id)
stan.data$id <- sapply(1:stan.data$N, function(i) which(ids == stan.data$id[i] ))

groups <- unique(stan.data$group)
stan.data$group <- sapply(1:stan.data$N, function(i) which(groups == stan.data$group[i] ))

stan.data$N_id <- length(unique(stan.data$id))
stan.data$N_groups <- length(unique(stan.data$group))


#save(stan.data, file = "standata211123")


library(rethinking)
#Baseline Multilevel Experience-weighted attraction model
m <- stan( file="Longitudinal_Conf.stan" , data=stan.data , chains=4, cores=4,refresh = 1, iter = 3000, control = list(adapt_delta=0.8, max_treedepth = 13))       
s <- extract.samples(m)


###
##
# Create plot (Fig. 4 in the manuscript)
##
###

#Crete color palette
par( mar = c(3,1.5,0,0.5), oma = c(1,2.5,2.5,0))

layout(matrix(c(1,1,1,1,
                2,2,   #gap
                3,3,
                4,4,4,4), 1, 12, byrow = TRUE))


keep[which(keep == -10)] <- 1
group[which(group == -10)] <- group[which(group == -10)+1]


col.pal <- colorRampPalette(brewer.pal(9, "Set1"))(N_groups)

#Time series of data
plot(1:30, xlim = c(1,30), ylim = c(1,stan.data$N_id), type = "n", ylab = "", xlab = "n")
for (i in 1:stan.data$N_id) {
  time <- max(which(keep[which(id == i)]==1))
  lines(1:time, rep(i, time), pch = 16, cex = 0.5, col = alpha("black", alpha=0.3))
  points(1:time, rep(i, time), pch = 16, cex = 0.5, col = alpha(col.pal[ group[which(id == i)][1:time] ], alpha=0.5))
}

mtext(side = 1, "Year", line = 2.5, cex = 1)
mtext(side = 2, "Participant ID", line = 2.5, cex = 1)
mtext('a', side=3, line=1, at=-1)

post <- inv_logit(s$logit_mu )
dens <- density(post)
x1 <- min(which(dens$x >= quantile(post, 0.05)))  
x2 <- max(which(dens$x <  quantile(post, 0.95)))
plot(dens, xlim = c(0.07,0.13), ylim = c(0,100), type="n", ann = FALSE, bty = "n", yaxt = "n")
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha(col.pal[1],alpha = 0.9), border = NA))

x1 <- min(which(dens$x >= quantile(post, 0)))  
x2 <- max(which(dens$x <  quantile(post, 1)))
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha(col.pal[1],alpha = 0.2), border = NA))

abline(v = mu, lty = 2, col = "black", lwd = 2)
mtext(side = 1, "Innovation rate", line = 2.5, cex = 1)
mtext('b', side=3, line=1, at=0.07)

post <- exp(s$log_theta)
dens <- density(post)
x1 <- min(which(dens$x >= quantile(post, 0.05)))  
x2 <- max(which(dens$x <  quantile(post, 0.95)))
plot(dens, xlim = c(2.5,4), ylim = c(0,5), type="n", ann = FALSE, bty = "n", yaxt = "n")
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha(col.pal[1],alpha = 0.9), border = NA))

x1 <- min(which(dens$x >= quantile(post, 0)))  
x2 <- max(which(dens$x <  quantile(post, 1)))
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha(col.pal[1],alpha = 0.2), border = NA))

abline(v = f, lty = 2, col = "black", lwd = 2)
mtext(side = 1, "Conformity exp.", line = 2.5, cex = 1)

MA <- stan.data$Max_age
plot(1:MA, type="n", ylim = c(0,0.8))

lower <- inv_logit( sapply(1:MA, function(i) HPDI(s$age_effects[,i],0.9))[1,])
upper <- inv_logit( sapply(1:MA, function(i) HPDI(s$age_effects[,i],0.9))[2,])

polygon(c(1:MA,MA:1), c(upper, rev(lower)), col=alpha(col.pal[1],alpha = 0.9), border = NA, ylim=c(0,5))

lower <- inv_logit( sapply(1:MA, function(i) HPDI(s$age_effects[,i],1))[1,])
upper <- inv_logit( sapply(1:MA, function(i) HPDI(s$age_effects[,i],1))[2,])

polygon(c(1:MA,MA:1), c(upper, rev(lower)), col=alpha(col.pal[1],alpha = 0.2), border = NA, ylim=c(0,5))

lines(age_mig_NL[2:MA], col = "black", lty = 2, lwd = 2)

mtext(side = 1, "Age", line = 2.5, cex = 1)
mtext(side = 2, "Migration rate", line = 2.5, cex = 1)
mtext('c', side=3, line=1, at=-5)

