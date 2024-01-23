library(parallel)

# initialize functions from ABM
source("functions_mig_conf.R")

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

png(filename = "abc_result.png", width = 18, height = 10, units = "cm", res = 500)

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

dev.off()

