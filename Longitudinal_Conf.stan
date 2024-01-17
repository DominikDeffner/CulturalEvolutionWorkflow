
//Multilevel time-series model to infer innovation rate and conformity
// D.Deffner, 2024, deffner@mpib-berlin.mpg.de

//Data block: Define and name the size of each observed variable
data{
   int N;              //Number of observations 
   int N_id;           //Number of individuals
   int N_groups;       //Number of groups
   int N_alt[N];           //Number of available traits for each choice
   int id[N];          //Unique individual identification
   int group[N];       //Group ID
   int choices[N];       //Chosen trait
   int innovate[N];       //Chosen trait
   real alternatives[N, 11]; //Matrix for all interaction partners
}

//Parameter block: Define and name the size of each unobserved variable.
parameters{
   real logit_mu;
   real log_theta;

   // Varying effects clustered on individual
    matrix[2,N_id] z_ID;
    vector<lower=0>[2] sigma_ID;
    cholesky_factor_corr[2] Rho_ID;

    // Varying effects clustered on groups
     matrix[2,N_groups] z_group;
     vector<lower=0>[2] sigma_group;
     cholesky_factor_corr[2] Rho_group;
}

//Transformed Parameters block: Here we multiply z-scores with variances and Cholesky factors to get varying effects back to right scale
transformed parameters{
      matrix[N_id,2] v_ID;
      matrix[N_groups,2] v_group;

      v_ID = ( diag_pre_multiply( sigma_ID , Rho_ID ) * z_ID )';
      v_group = ( diag_pre_multiply( sigma_group , Rho_group ) * z_group )';
}

//Model block: Here compute the log posterior
model{

  //Priors
   logit_mu  ~ normal(0,1);
   log_theta ~ normal(0,1);

  //Varying effects priors
  to_vector(z_ID) ~ normal(0,1);
  sigma_ID ~ exponential(1);
  Rho_ID ~ lkj_corr_cholesky(4);

  to_vector(z_group) ~ normal(0,1);
  sigma_group ~ exponential(1);
  Rho_group ~ lkj_corr_cholesky(4);

//For each choice we first estimate the probability that an individual innovates
//If they didn't innovate, we also estimate the strength of (anti)conformity
for (i in 1:N){

//Probability of innovation
target += bernoulli_logit_lpmf(innovate[i]| logit_mu + v_ID[id[i], 1] + v_group[group[i], 1] );

if (innovate[i] == 0){
  //Vector for choice probabilities
  vector[N_alt[i]] p;

  //Compute choice probabilities based on individual- and group-specific conformity value
   for ( j in 1:N_alt[i] ) p[j] = alternatives[i,j]^exp(log_theta + v_ID[id[i], 2] + v_group[group[i], 2]);
    p = p / sum(p);

  //Add log probability of observed trait choice to target
  target += categorical_lpmf(choices[i] | softmax(p));
  
  }

}

}// end model