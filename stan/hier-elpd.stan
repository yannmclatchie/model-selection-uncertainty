//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//

data {
  // Data dimensions
  int<lower=1> N; // Number of datapoints
  int<lower=2> K; // Number of models
  int<lower=1> I; // Number of observations
  
  vector[N] y; // Point-wise elpd values stacked into a vector
  array[N] int<lower=1, upper=K> models; // Model indicies
  array[N] int<lower=1, upper=I> obs; // Observation indicies
}

transformed data {
  real mean_y;
  real<lower=0> sd_y;
  
  mean_y = mean(y);
  sd_y = sd(y);
}

parameters {
  real alpha;
  real<lower=0> sigma;
  
  // Model level effects
  vector[K] beta_raw;
  real mu_beta;
  real<lower=0> sigma_beta;
  
  // Datapoint level effects
  vector[I] gamma_raw;
  real mu_gamma;
  real<lower=0> sigma_gamma;
  
  // Interactions
  matrix[K, I] delta_raw;
  real mu_delta;
  real<lower=0> sigma_delta;
}

transformed parameters {
  vector[N] mu;
  vector[K] beta;
  vector[I] gamma;
  matrix[K, I] delta;
  
  beta = mu_beta + beta_raw * sigma_beta;
  gamma = mu_gamma + gamma_raw * sigma_gamma;
  
  delta = rep_matrix(mu_delta, K, I) + delta_raw * sigma_delta;

  for(n in 1:N){
    mu[n] = alpha + beta[models[n]] + gamma[obs[n]] + delta[models[n], obs[n]];
  }
}

model {
  // Intercept
  alpha ~ normal(mean_y, sd_y * 2.5);

  // Beta
  sigma_beta ~ normal(0, 0.5);
  mu_beta ~ normal(0, 0.5);
  beta_raw ~ std_normal();
  
  // Gamma
  sigma_gamma ~ cauchy(0, 0.5);
  mu_gamma ~ normal(0, 0.5);
  gamma_raw ~ std_normal();
  
  // Delta
  sigma_delta ~ cauchy(0, 0.5);
  mu_delta ~ normal(0, 0.5);
  for(k in 1:K){
    for(n in 1:I){
      delta_raw[k, n] ~ std_normal();   
    }
  }
  
  // Likelihood
  sigma ~ normal(0, 1);
  y ~ normal(mu, sigma); 
} 

generated quantities {
  vector[N] y_rep;
  for(n in 1:N){
    y_rep[n] = normal_rng(
      alpha + beta[models[n]] + gamma[obs[n]] + delta[models[n], obs[n]],
      sigma);
  }
}
