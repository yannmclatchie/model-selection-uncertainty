data {
  // dimensions
  int<lower=0> N_train;
  int<lower=0> N_test;
  int<lower=0> d;
  // predictor matrices
  matrix[N_train, d] x_train;
  matrix[N_test, d] x_test;
  // response vectors
  vector[N_train] y_train;
  vector[N_test] y_test;
}
parameters {
  vector[d] beta;
  real<lower=0> sigma;
}
model {
  // priors
  beta ~ normal(0, 1);
  sigma ~ normal(0, 1);
  // likelihood
  y_train ~ normal(x_train * beta, sigma);
}
generated quantities {
  vector[N_train] log_lik;
  vector[N_test] log_lik_test;
  // in-sample log-likelihood
  for (n in 1:N_train) {
    log_lik[n] = normal_lpdf(y_train[n] | x_train[n] * beta, sigma);
  }
  // test log-likelihood
  for (n in 1:N_test) {
    log_lik_test[n] = normal_lpdf(y_test[n] | x_test[n] * beta, sigma);
  }
}
