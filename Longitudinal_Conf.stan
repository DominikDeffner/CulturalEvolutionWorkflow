
// Deffner, D., Fedorova, N., Andrews, J. & McElreath, R.
// Bridging theory and data: A computational workflow for cultural evolution
// Longitudinal transmission analysis to infer conformity, innovation and age-dependent migration (authored by D.Deffner; email: deffner@mpib-berlin.mpg.de)


functions{
//Function for Gaussian Process kernel
  matrix GPL(int K, real C, real D, real S){
   matrix[K,K] Rho;
   real KR;
   KR = K;
   for(i in 1:(K-1)){
   for(j in (i+1):K){
    Rho[i,j] = C * exp(-D * ( (j-i)^2 / KR^2) );
    Rho[j,i] = Rho[i,j];
    }}
   for (i in 1:K){
    Rho[i,i] = 1;
    }
   return S*cholesky_decompose(Rho);
  }

}

//Data block: Define and name the size of each observed variable
data{
   int N;              //Number of observations 
   int N_id;           //Number of individuals
   int N_groups;       //Number of groups
   int N_partners;     //Number of interaction partners
   int N_alt[N];       //Number of available traits for each choice
   int id[N];          //Unique individual identification
   int group[N];       //Group ID
   int age[N];         //Age
   int Max_age;        //maximum age
   int choices[N];     //Chosen trait
   int innovate[N];    //Does individual have a new variant?
   int migrate[N];    //Did individual migrate since the last timestep?
   real frequencies[N, N_partners+1]; //Matrix for all interaction partners
}

//Parameter block: Define and name the size of each unobserved variable.
parameters{
  #Innovation rate, average migration rate and conformity on the latent (logit/log) scale
   real logit_mu;
   real logit_m;
   real log_theta;
   
   //Vector for Gaussian process age effects
   vector[Max_age] age_offsets;    

  //Here we define the Control parameters for the Gaussian processes; they determine how covariance changes with increasing distance in age
  real<lower=0> eta;
  real<lower=0> sigma;
  real<lower=0, upper=1> rho;

   // Varying effects clustered on individual
    matrix[3,N_id] z_ID;
    vector<lower=0>[3] sigma_ID;
    cholesky_factor_corr[3] Rho_ID;

    // Varying effects clustered on groups
     matrix[3,N_groups] z_group;
     vector<lower=0>[3] sigma_group;
     cholesky_factor_corr[3] Rho_group;
}

//Transformed Parameters block: Here we multiply z-scores with variances and Cholesky factors to get varying effects back to right scale
transformed parameters{
      matrix[N_id,3] v_ID;
      matrix[N_groups,3] v_group;
      vector[Max_age] age_effects; 
      v_ID = ( diag_pre_multiply( sigma_ID , Rho_ID ) * z_ID )';
      v_group = ( diag_pre_multiply( sigma_group , Rho_group ) * z_group )';
      
      age_effects = logit_m + age_offsets;
}

//Model block: Here compute the log posterior
model{

  //Priors
   logit_mu  ~ normal(0,2);
   logit_m  ~ normal(0,2);
   log_theta ~ normal(0,1);
   
   //Gaussian process on migration rate
  eta ~ exponential(100);
  sigma ~ exponential(1);
  rho ~ beta(300, 1);

  //Varying effects priors
  to_vector(z_ID) ~ normal(0,1);
  sigma_ID ~ exponential(1);
  Rho_ID ~ lkj_corr_cholesky(4);

  to_vector(z_group) ~ normal(0,1);
  sigma_group ~ exponential(1);
  Rho_group ~ lkj_corr_cholesky(4);

  //We compute age-specific offsets
  age_offsets ~ multi_normal_cholesky( rep_vector(0, Max_age) , GPL(Max_age, rho, eta, sigma) );
  
//For each choice we first estimate the probability that an individual innovates (and migrates)
//If they didn't innovate, we also estimate the strength of (anti)conformity
for (i in 1:N){

//Probability of innovation
target += bernoulli_logit_lpmf(innovate[i]| logit_mu + v_ID[id[i], 1] + v_group[group[i], 1] );

//Probability of migration
target += bernoulli_logit_lpmf(migrate[i]| age_effects[age[i]] + v_ID[id[i], 2] + v_group[group[i], 2] );

if (innovate[i] == 0){
  //Vector for choice probabilities
  vector[N_alt[i]] p;

  //Compute choice probabilities based on individual- and group-specific conformity value
   for ( j in 1:N_alt[i] ) p[j] = frequencies[i,j]^exp(log_theta + v_ID[id[i], 3] + v_group[group[i], 3]);
    p = p / sum(p);

  //Add log probability of observed trait choice to target
  target += categorical_lpmf(choices[i] | p);
  }

}

}// end model
