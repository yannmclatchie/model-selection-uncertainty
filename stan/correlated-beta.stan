//
// This Stan program defines a the correlated effect
// model over model elpds.
//

data {
  int<lower=1> N;
  int<lower=2> K;
  vector[K] y[N];
}

parameters {
  // base elpd level
  real common_elpd;
  
  // model level effects
  vector[K] alpha;
  vector[K] mu_beta;
  vector<lower=0>[K] sigma_beta;
  cholesky_factor_corr[K] Lcorr;  
  
  // residual
  real<lower=0> sigma;
}

transformed parameters {
  // model effects
  vector[K] beta;
  beta = mu_beta + sigma_beta .* (Lcorr * alpha);
}

model {
  // priors
  common_elpd ~ std_normal();
  sigma ~ cauchy(0, 5); 
  alpha ~ std_normal();
  sigma_beta ~ cauchy(0, 5); 
  mu_beta ~ std_normal();

  // likelihood
  for (n in 1:N) {
    for (k in 1:K) {
      y[n,k] ~ normal(common_elpd + beta[k], sigma);
    }
  }
}

generated quantities {
  // compute the covariance and correlation matrices 
  matrix[K,K] Omega;
  matrix[K,K] Sigma;
  Omega = multiply_lower_tri_self_transpose(Lcorr);
  Sigma = quad_form_diag(Omega, sigma_beta); 
  
  // likelihood
  vector[K] log_lik[N];
  for (n in 1:N) {
    for (k in 1:K) {
      log_lik[n,k] = normal_lpdf(y[n,k] | common_elpd + beta[k], sigma);
    }
  }
  
  // joint prior specification
  real lprior;
  lprior = normal_lpdf(common_elpd | 0, 1) +
    cauchy_lpdf(sigma | 0, 5) +
    normal_lpdf(alpha | 0, 1) +
    cauchy_lpdf(sigma_beta | 0, 5) +
    normal_lpdf(mu_beta | 0, 1);
}
