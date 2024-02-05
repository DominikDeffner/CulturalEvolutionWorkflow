# Cultural Evolution Workflow

This repository contains the scripts to reproduce all results and plots for the workflow examples in 

***Deffner, D., Fedorova, N., Andrews, J. & McElreath, R. (submitted) Bridging theory and data: A computational workflow for cultural evolution***

***Workflow examples***
- Moving from DAGS to data: "DAGstoData.R" simulates synthetic data for 30 different countries based on the DAG in Fig.2b and performs a power analysis varying the amount of unexplained variation and the effect size.
- Agent-based simulations: "ConDiv_ABM.R" simulates an agent-based model combining migration and conformity and plots the results (Fig.3 in the manuscript).
- Longitudinal transmission analysis: "ConDiv_LongitudinalTransmission.R" simulates longitudinal data on cultural traits and social networks for hypothetical participants. It then sources "Longitudinal_Conf.stan" and fits a time-series transmission model in stan and plots the results (Fig.4 in the manuscript). The model uses real-world age trajectories for migration stored in "beta_df.RDS". 
- Approximate Bayesian Computation (ABC): "ABC_analysis.R" simulates cross-sectional data on cultural Fst for hypothetical participants. It then uses Approximate Bayesian Computation to infer parameters and plots the results (Fig.5 in the manuscript)

***Software requirements***
The code was written in R 4.0.3. Statistical models are fit using the Stan MCMC engine via the rstan package (2.21.2), which requires a C++ compiler. Installation instructions are available at https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started. See also the Stan user guide at https://mc-stan.org/users/documentation. The rethinking package (2.12) is required to process fitted model outputs (installation instructions at http://xcelab.net/rm/software/).


