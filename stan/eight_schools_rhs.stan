functions {
  /* Efficient computation of the horseshoe prior
   * see Appendix C.1 in https://projecteuclid.org/euclid.ejs/1513306866
   * Args:
   *   z: standardized population-level coefficients
   *   lambda: local shrinkage parameters
   *   tau: global shrinkage parameter
   *   c2: slab regularization parameter
   * Returns:
   *   population-level coefficients following the horseshoe prior
   */
  vector horseshoe(vector z, vector lambda, real tau, real c2) {
    int K = rows(z);
    vector[K] lambda2 = square(lambda);
    vector[K] lambda_tilde = sqrt(c2 * lambda2 ./ (c2 + tau ^ 2 * lambda2));
    return z .* lambda_tilde * tau;
  }
}
data {
  int<lower=0> J; // number of schools
  array[J] real y; // estimated treatment effect (school j)
  array[J] real<lower=0> sigma; // std err of effect estimate (school j)
  
  // data for the horseshoe prior
  real<lower=0> hs_df; // local degrees of freedom
  real<lower=0> hs_df_global; // global degrees of freedom
  real<lower=0> hs_df_slab; // slab degrees of freedom
  real<lower=0> hs_scale_global; // global prior scale
  real<lower=0> hs_scale_slab; // slab prior scale
}
parameters {
  real mu;
  
  // local parameters for the horseshoe prior
  vector[J] phi_z;
  vector<lower=0>[J] hs_local;
  // horseshoe shrinkage parameters
  real<lower=0> hs_global; // global shrinkage parameter
  real<lower=0> hs_slab; // slab regularization parameter
  real<lower=0> lambda; // dispersion parameter
}
transformed parameters {
  // regularised horseshoe likelihood
  vector[J] theta = horseshoe(phi_z, hs_local, hs_global, hs_scale_slab ^ 2 * hs_slab);
}
model {
  // horseshoe priors
  lambda ~ std_normal();
  phi_z ~ std_normal();
  hs_global ~ student_t(hs_df_global, 0, hs_scale_global * lambda);
  hs_slab ~ inv_gamma(0.5 * hs_df_slab, 0.5 * hs_df_slab);
  hs_local ~ student_t(hs_df, 0, 1);
  
  // likelihood
  y ~ normal(theta, sigma);
}
