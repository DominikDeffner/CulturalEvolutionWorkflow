
###########################################################################################
## Deffner, D., Fedorova, N., Andrews, J. & McElreath, R.
## Bridging theory and data: A computational workflow for cultural evolution
## Approximate Bayesian Computation (authored by N.Fedorova; email: natalia_fedorova@eva.mpg.de)
###########################################################################################

#This script simulates cross-sectional data on cultural Fst for hypothetical participants
#It then uses Approximate Bayesian Computation to infer parameters and plots the results (Fig.5 in the manuscript)

#Set working directory to GitHub repo
setwd("~/GitHub/CulturalEvolutionWorkflow")

library(parallel)

# initialize functions from ABM


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



# generate reference Fst - i.e. data we collected in the field
# here we generate reference data that comes from a population where transmission is unbiased
reference_data_mig_un <- mig_abm(N = 3000,         
                              N_groups = N_groups,      
                              N_steps = 100,    
                              N_burn_in = 1000,
                              N_mod = 30,       
                              m_NL = FALSE,    
                              m_const = 0.3,
                              mu = 0.1,       
                              theta = 1,
                              r_dist  = 0.2)


# vary conformity(theta) and migration rate(m_const) - we assume flat priors
# create jobs with param combinations to be entered into mig_abm

# number of combinations needed for ABC analysis
n_comb <- 10

# params to input
N_un <- rep(3000, n_comb)
N_groups_un <- rep(N_groups, n_comb)  
N_steps_un <- rep(100, n_comb) 
N_burn_in_un <- rep(1000, n_comb)
N_mod_un <- rep(30, n_comb)
m_NL_un <- rep(FALSE, n_comb)

m_seq_un <- seq(0, 1, 0.1)
m_const_un <- sample(m_seq_un, n_comb, replace = TRUE)

mu_un <- rep(0.05, n_comb)

theta_seq_un <- c(0.5, 0, 1, 3)
theta_un <- sample(theta_seq_un, n_comb, replace = TRUE)

r_dist_un <- rep(0.2, n_comb)
fst_compare_un <- rep(TRUE, n_comb)

jobs_fst_compare_mig_un <- data.frame(N_un, 
                                   N_groups_un, 
                                   N_steps_un, 
                                   N_burn_in_un,
                                   N_mod_un,
                                   m_NL_un,
                                   m_const_un,
                                   mu_un,
                                   theta_un,
                                   r_dist_un,
                                   fst_compare_un)


cores <- 80

# ABC loop which compares Fst values from reference df and sim df
abc_fst_diff_mig_un <- mclapply(1:n_comb,
                             function(i) {
                               if (i %% 100 == 0) print(i)
                               out <- abc_loop_fst(
                                 N = jobs_fst_compare_mig_un$N_un[i],         
                                 N_groups = jobs_fst_compare_mig_un$N_groups_un[i],      
                                 N_steps = jobs_fst_compare_mig_un$N_steps_un[i],    
                                 N_burn_in = jobs_fst_compare_mig_un$N_burn_in_un[i],
                                 N_mod = jobs_fst_compare_mig_un$N_mod_un[i],       
                                 m_NL = jobs_fst_compare_mig_un$m_NL_un[i],    
                                 m_const = jobs_fst_compare_mig_un$m_const_un[i],
                                 mu = jobs_fst_compare_mig_un$mu_un[i],       
                                 theta = jobs_fst_compare_mig_un$theta_un[i], 
                                 r_dist  = jobs_fst_compare_mig_un$r_dist_un[i],
                                 reference_data = reference_data_mig_un,
                                 fst_compare = jobs_fst_compare_mig_un$fst_compare_un[i])},
                             mc.cores = cores)


# save abc_output with jobs and reference_data
abc_output_fst_mig_un <- list(abc_fst_diff_mig_un, jobs_fst_compare_mig_un, reference_data_mig_un)
saveRDS(abc_output_fst_mig_un, file="abc_output_fst_mig_un.RDS")


### ABC analysis for conformist transmission
# generate reference Fst - i.e. data we collected in the field
# here we generate reference data that comes from a population where transmission is conformist
reference_data_mig_c <- mig_abm(N = 3000,         
                                 N_groups = N_groups,      
                                 N_steps = 100,    
                                 N_burn_in = 1000,
                                 N_mod = 30,       
                                 m_NL = FALSE,    
                                 m_const = 0.5,
                                 mu = 0.1,       
                                 theta = 3,
                                 r_dist  = 0.2)


# vary conformity(theta) and migration rate(m_const) - we assume flat priors
# create jobs with param combinations to be entered into mig_abm
# note: n_comb and cores defined above

# params to input
N_c <- rep(3000, n_comb)
N_groups_c <- rep(N_groups, n_comb)  
N_steps_c <- rep(100, n_comb) 
N_burn_in_c <- rep(1000, n_comb)
N_mod_c <- rep(30, n_comb)
m_NL_c <- rep(FALSE, n_comb)

m_seq_c <- seq(0, 1, 0.1)
m_const_c <- sample(m_seq_c, n_comb, replace = TRUE)

mu_c <- rep(0.05, n_comb)

theta_seq_c <- c(0.5, 0, 1, 3)
theta_c <- sample(theta_seq_c, n_comb, replace = TRUE)

r_dist_c <- rep(0.2, n_comb)
fst_compare_c <- rep(TRUE, n_comb)

jobs_fst_compare_mig_c <- data.frame(N_c, 
                                   N_groups_c, 
                                   N_steps_c, 
                                   N_burn_in_c,
                                   N_mod_c,
                                   m_NL_c,
                                   m_const_c,
                                   mu_c,
                                   theta_c,
                                   r_dist_c,
                                   fst_compare_c)


# ABC loop which compares Fst values from reference df and sim df
abc_fst_diff_mig_c <- mclapply(1:n_comb,
                             function(i) {
                               if (i %% 100 == 0) print(i)
                               out <- abc_loop_fst(
                                 N = jobs_fst_compare_mig_c$N_c[i],         
                                 N_groups = jobs_fst_compare_mig_c$N_groups_c[i],      
                                 N_steps = jobs_fst_compare_mig_c$N_steps_c[i],    
                                 N_burn_in = jobs_fst_compare_mig_c$N_burn_in_c[i],
                                 N_mod = jobs_fst_compare_mig_c$N_mod_c[i],       
                                 m_NL = jobs_fst_compare_mig_c$m_NL_c[i],    
                                 m_const = jobs_fst_compare_mig_c$m_const_c[i],
                                 mu = jobs_fst_compare_mig_c$mu_c[i],       
                                 theta = jobs_fst_compare_mig_c$theta_c[i], 
                                 r_dist  = jobs_fst_compare_mig_c$r_dist_c[i],
                                 reference_data = reference_data_mig_c,
                                 fst_compare = jobs_fst_compare_mig_c$fst_compare_c[i])},
                             mc.cores = cores)


# save abc_output with jobs and reference_data
abc_output_fst_mig_c <- list(abc_fst_diff_mig_c, jobs_fst_compare_mig_c, reference_data_mig_c)
saveRDS(abc_output_fst_mig_c, file="abc_output_fst_mig_c.RDS")


### Extracting posteriors with rejection algorithm

# if loading abc result, uncomment here:
#abc_output_fst_mig_un <- readRDS("abc_output_fst_mig.RDS")
#reference_data_un <- abc_output_fst_mig[[3]]

# Otherwise:
# rejection algorithm
# sort parameter combinations by fst differences, select 1000 "best" parameter combinations
# this constitutes the joint posterior

#extract jobs - parameter combinations used to generate data in ABC
comb_test_un <- as.data.frame(abc_output_fst_mig_un[[2]]) 
#extract fst differences from ABC algorithm, and make them absolute
comb_test_un$abc_diff <- unlist(abc_output_fst_mig_un[[1]])
comb_test_un$absolute_diff <- abs(comb_test_un$abc_diff)
#order by differences
output_ordered_abs_un <- comb_test_un[order(comb_test_un$absolute_diff, decreasing = FALSE),]

#extract best 1000
abc_posterior_fst_mig_un <- output_ordered_abs_un[1:1000,]


# and then the posterior for reference data with conformist transmission

# if loading abc result, uncomment lines:
#abc_output_fst_mig_c <- readRDS("abc_output_fst_mig_c.RDS")
#reference_data_c <- abc_output_fst_mig_2[[3]]

#Otherwise
comb_test_c <- as.data.frame(abc_output_fst_mig_c[[2]]) 
comb_test_c$abc_diff <- unlist(abc_output_fst_mig_c[[1]])
comb_test_c$absolute_diff <- abs(comb_test_c$abc_diff)
output_ordered_abs_c <- comb_test_c[order(comb_test_c$absolute_diff, decreasing = FALSE),]
abc_posterior_fst_mig_c <- output_ordered_abs_c[1:1000,]


### posterior prediction
# to view what the joint posterior implies on the outcome scale, here in terms of fst, we need to push posterior
# estimates back through the model

# posterior prediction for reference data with unbiased transmission

# posterior predictions to generate:
post_pred_n <- 10
post_pred_un <- list()

for(i in 1:post_pred_n){
  post_pred_un[[i]] <- mig_abm(   N = abc_posterior_fst_mig_un$N_un[i],         
                                  N_groups = abc_posterior_fst_mig_un$N_groups_un[i],      
                                  N_steps = abc_posterior_fst_mig_un$N_steps_un[i],    
                                  N_burn_in = abc_posterior_fst_mig_un$N_burn_in_un[i],
                                  N_mod = abc_posterior_fst_mig_un$N_mod_un[i],       
                                  m_NL = abc_posterior_fst_mig_un$m_NL_un[i],    
                                  m_const = abc_posterior_fst_mig_un$m_const_un[i],
                                  mu = abc_posterior_fst_mig_un$mu_un[i],       
                                  theta = abc_posterior_fst_mig_un$theta_un[i], 
                                  r_dist  = abc_posterior_fst_mig_un$r_dist_un[i])
  
}

# select fst difference for the 100th timestep from all runs
last_fst_un <- unlist(lapply(post_pred_un, function(l) l[[1]][100]))



# posterior prediction for reference data with conformist transmission
# note: number to generate defined above

post_pred_c <- list()

for(i in 1:post_pred_n){
  post_pred_c[[i]] <- mig_abm(   N = abc_posterior_fst_mig_c$N_c[i],         
                                  N_groups = abc_posterior_fst_mig_c$N_groups_c[i],      
                                  N_steps = abc_posterior_fst_mig_c$N_steps_c[i],    
                                  N_burn_in = abc_posterior_fst_mig_c$N_burn_in_c[i],
                                  N_mod = abc_posterior_fst_mig_c$N_mod_c[i],       
                                  m_NL = abc_posterior_fst_mig_c$m_NL_c[i],    
                                  m_const = abc_posterior_fst_mig_c$m_const_c[i],
                                  mu = abc_posterior_fst_mig_c$mu_c[i],       
                                  theta = abc_posterior_fst_mig_c$theta_c[i], 
                                  r_dist  = abc_posterior_fst_mig_c$r_dist_c[i])
  
}

# select fst difference for the 100th timestep from all runs
last_fst_c = unlist(lapply(post_pred_c, function(l) l[[1]][100]))

### plotting results

# extract theta posteriors for both runs for plotting
# note: plotting has been tailored to default runs described in this script
# xlim and ylim as well as other plotting parameters may need altering for different runs

theta_post_un <- as.data.frame(table(abc_posterior_fst_mig_un$theta_un), stringsAsFactors = FALSE)
theta_post_un <- rbind(theta_post_un, c(2,0))
theta_post_un <- rbind(theta_post_un, c(3,0))

theta_post_c <- as.data.frame(table(abc_posterior_fst_mig_c$theta_c), stringsAsFactors = FALSE)
theta_post_c <- rbind(theta_post_c, c(0,0))
theta_post_c <- rbind(theta_post_c, c(0.5,0))
theta_post_c <- rbind(theta_post_c, c(1,0))
theta_post_c <- theta_post_c[order(theta_post_c$Var1),]

# plotting parameters
ylim_i <- 600
at_l <- -0.5

#png(filename = "abc_result.png", width = 18, height = 10, units = "cm", res = 500)

par(mfrow = c(2,3),
    mar = c(4,4,3,3),
    oma = c(2,2,1,2))

plot(theta_post_un,
     type = "b",
     #main = "Conformity exp. posterior", 
     lwd = 2,
     xlab = "Conformity exp.",
     ylab = "Frequency",
     ylim = c(0,ylim_i),
     xlim = c(0,3),
     col = "red3",
     bty = "n")
abline(v = 1, lty = 3)
mtext("a", side = 3, line = 1, at = at_l)

m_const_post_un <- as.data.frame(table(abc_posterior_fst_mig_un$m_const_un), stringsAsFactors = FALSE)
plot(m_const_post_un,
     type = "b",
     #main = "Migration rate posterior", 
     lwd = 2,
     xlab = "Migration rate",
     ylab = "Frequency",
     ylim = c(0,ylim_i),
     xlim = c(0,0.5),
     col = "red3",
     bty = "n")
abline(v = 0.3, lty = 3)
mtext("b", side = 3, line = 1, at = -0.1)


plot(density(last_fst_un),
     bty = "n",
     col = "red3",
     lwd = 2,
     main = "",
     xlim = c(0.019,0.025),
     xlab = "Cultural Fst")
abline(v = reference_data_mig_un[[1]][100], lty = 2, lwd = 2)
mtext("c", side = 3, line = 1, at = 0.0175)

### different ref

plot(theta_post_c,
     type = "b",
     #main = "Conformity exp. posterior", 
     lwd = 2,
     xlab = "Conformity exp.",
     ylab = "Frequency",
     ylim = c(0,800),
     xlim = c(0,3),
     col = "red3",
     bty = "n")
abline(v = 3, lty = 3)
mtext("d", side = 3, line = 1, at = at_l)

m_const_post_c <- as.data.frame(table(abc_posterior_fst_mig_c$m_const_c), stringsAsFactors = FALSE)
plot(m_const_post_c,
     type = "b",
     #main = "Migration rate posterior", 
     lwd = 2,
     xlab = "Migration rate",
     ylab = "Frequency",
     ylim = c(0,200),
     xlim = c(0,0.6),
     col = "red3",
     bty = "n")
abline(v = 0.5, lty = 3)
mtext("e", side = 3, line = 1, at = -0.1)


plot(density(last_fst_c),
     bty = "n",
     col = "red3",
     lwd = 2,
     main = "",
     xlab = "Cultural Fst")
abline(v = reference_data_mig_c[[1]][100], lty = 2, lwd = 2)
mtext("f", side = 3, line = 1, at = 0.28)

#dev.off()

