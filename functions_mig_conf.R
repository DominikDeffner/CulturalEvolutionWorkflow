
# Cultural evolution model with age-dependent learning and migration
# We assume a meta-population of different groups (or villages) which are connected through migration.
# Individuals can learn cultural traits throughout their lives but become less likely to update
# their traits as they get older


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


#Spatial components
#Grid size for spatially-explicit model with local migration
size_grid <- 50

#Number of groups
N_groups <- 30

#First create rectangular grid for villages, we try to make it spatially explicit
grid <- matrix(0, nrow = size_grid, ncol = size_grid)

# Now let's populate the grid with villages
Locations <- sample(size_grid * size_grid, size = N_groups)
grid[Locations] <- sample(N_groups)



#'mig_abm 
#' 
# parameters:
# N = 3000             #Number of agents
# N_groups = N_groups         #Number of groups
# N_steps = 100    #Number of time steps
# N_burn_in = 1000 #Number of time steps to get equilibrium age structure
# N_mod = 30       #Number of role models
# m_NL = TRUE     #If TRUE, we take real migration rates from NL; if FALSE constant rate m_const
# m_const = 0.1
# mu = 0.1        #Innovation rate (fraction of learning events that are innovations)
# theta = 1.3            #Conformity parameter
# dist = dist        # distance matrix for villages
# r_dist  = 0.2      #Effect of distance on migration (if r_dist = 0, migrants choose village randomly, 0.1 seems to be good value for local migration)
#'
#' @return Fst = a vector of Fst values for the population at each timestep of N_steps

mig_abm <- function(N = 3000,         
                    N_groups = N_groups,      
                    N_steps = 100,    
                    N_burn_in = 1000,
                    N_mod = 30,       
                    m_NL = TRUE,    
                    m_const = 0.1,
                    mu = 0.1,       
                    theta = 1.3,
                    r_dist  = 0.2   
){
  
  # number of inds per group
  N_per_group <- N/N_groups
  
  #Define developmental rates for learning and mortality
  r_learn = 0.03
  r_mort  = 0.001
  
  #Assign migration rates
  if (m_NL == TRUE){
    m <- age_mig_NL
  } else {
    # m <- rep( sum(age_mig_NL*Age_hat), max_age )
    m <- rep(m_const, max_age)
  }
  
  # Initialize population; we just keep track of personal id, group id and conformity exponent per individual
  Age <- sample(1:80, N, replace = TRUE)
  group <- rep(1:N_groups,each= N_per_group)
  
  # calculate distance between villages
  dist <- matrix(0, nrow = N_groups, ncol = N_groups)
  for (x in 1:N_groups) {
    for (y in 1:N_groups) {
      dist[x,y] <- sqrt( (which(grid == x, arr.ind = TRUE)[1] - which(grid == y, arr.ind = TRUE)[1])^2 +
                           (which(grid == x, arr.ind = TRUE)[2] - which(grid == y, arr.ind = TRUE)[2])^2 )
    }
  }
  
  #Initialize cultural traits 
  Traits <- rep(0, N)
  
  #We start with unique variants per group
  unique_per_group <- sample(1:N_groups,  replace = FALSE)
  for (g in 1:N_groups) Traits[group == g] <- unique_per_group[g]
  
  #Counter for variants to make sure that innovations are really new
  Counter <- max(unique_per_group)
  
  #Initialize Fst vector
  Fst <- rep(0, N_steps)
  total.var <- rep(0, N_steps)
  between.var <- rep(0, N_steps)
  
  
  #Loop over years
  for (t in 1: (N_burn_in + N_steps ) ){
    #print(t)
    
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
      
      #Loop over all learners
      
      for (i in Learners) {
        
        #Sample models
        Model_ids <- sample(which(group == group[i]), N_mod)
        
        #Vector with unique variants
        Variants <- unique(Traits[c(Model_ids, i)]) 
        
        #Frequency of each variant
        Freq_Variants <- c()
        for (x in Variants) Freq_Variants[which(Variants == x)] <- length(which(Traits[Model_ids] == x))
        
        #Probability individuals choose each variant
        prob <- Freq_Variants^theta / sum(Freq_Variants^theta)
        
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
      }#i
      
      #Replace old by new cultural traits
      Traits[Learners] <- Traits_new
      
      #Function for getting fixation Index (Fst) from trait vectors (based on Mesoudi, 2018, Migration, acculturation, and the maintenance of between-group cultural variation)
      J <- unique(Traits)
      frequencies <- matrix(0, nrow = N_groups, ncol = length(J))
      for (g in 1:N_groups) frequencies[g,] <- sapply(J, function(j) length(which(Traits[which(group==g)] == j)) / N_per_group )
      
      total.var_t <- 1 - sum(colMeans(frequencies)^2)  # 1 - sum of squared means of each trait
      within.var_t <- mean(1 - rowSums(frequencies^2)) # mean of each group's (1 - the sum of squared freq of each trait)
      between.var_t <- total.var_t - within.var_t
      fst_t <- between.var_t / total.var_t  #Fst
      
      #Quantify Fst
      Fst[t- N_burn_in] <- fst_t
      total.var[t- N_burn_in] <- total.var_t
      between.var[t- N_burn_in] <- between.var_t
      
    }
    
  }#t
  
  return(list(Fst, total.var, between.var))
  
}#mig_abm function




#' abc_loop_fst
#' 
# N = 3000                    #Number of agents in ABM
# N_groups = N_groups         #Number of groups in ABM
# N_steps = 100               #Number of time steps
# N_burn_in = 1000            #Number of time steps to get equilibrium age structure
# N_mod = 30                  #Number of role models
# m_NL = TRUE                 #If TRUE, we take real migration rates from NL; if FALSE constant rate m_const
# m_const = 0.1
# mu = 0.1                    #Innovation rate (fraction of learning events that are innovations)
# theta = 1.3                 #Conformity parameter theta
# dist = dist                 #distance matrix for villages
# r_dist  = 0.2               #Effect of distance on migration (if r_dist = 0, migrants choose village randomly, 0.1 seems to be good value for local migration)
#' reference_data             #data that is being studied
#' fst_compare = TRUE         # TRUE = Fst values compred, FALSE = compare within and between group variance
#'
#' @return fst_diff

abc_loop_fst <- function(N = 3000,         
                         N_groups = N_groups,      
                         N_steps = 100,    
                         N_burn_in = 1000,
                         N_mod = 30,       
                         m_NL = TRUE,    
                         m_const = 0.1,
                         mu = seq(0, 0.1, 0.01),       
                         theta = seq(0.5, 5, 0.5),           
                         r_dist  = 0.2,
                         reference_data = reference_data,
                         fst_compare = TRUE){
  
  # generate data from abm to compare
  sim_data <- mig_abm(N = N,         
                      N_groups = N_groups,      
                      N_steps = N_steps,    
                      N_burn_in = N_burn_in,
                      N_mod = N_mod,       
                      m_NL = m_NL,    
                      m_const = m_const,
                      mu = mu,       
                      theta = theta,           
                      r_dist  = r_dist)
  
  print(fst_compare)
  print(mu)
  print(theta)
  print(sim_data[[1]][N_steps])
  print(reference_data[[1]][N_steps])
  
  if(fst_compare == TRUE){
    # compare reference and sim data on fst 
    # calculate Fst difference between sim and reference data
    fst_diff <- sim_data[[1]][N_steps] - reference_data[[1]][N_steps]
    output <- fst_diff
    print(fst_diff)
    
  }else{
    # compare reference and sim data based on total and between var separately
    total.var_diff <- sim_data[[2]][N_steps] - reference_data[[2]][N_steps]
    between.var_diff <- sim_data[[3]][N_steps] - reference_data[[3]][N_steps]
    
    output <- abs(total.var_diff) + abs(between.var_diff)
    
  }
  
  return(output)
  
}


