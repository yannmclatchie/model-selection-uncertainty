library(cmdstanr)
library(simstudy)
library(loo)
library(ggplot2)
library(geomtextpath)
library(dplyr)
library(purrr)
stanfit <- function(fit) rstan::read_stan_csv(fit$output_files())

simulate_data <- function(n, K, eps, beta_delta) {
  # define the DGP 
  def <- defData(varname = "x0", formula = "1")
  def <- defRepeat(def, nVars = K, prefix = "x", formula = "0",
                   variance = "1", dist = "normal")
  def <- defData(def, "y", formula = "0 * x0 + 1 + ..beta_delta * x1", 
                   variance = "..eps", dist = "normal")
  
  # generate the data
  dd_train <- genData(n, def)
  dd_test <- genData(n, def)
  
  # return output
  return(list(train = dd_train, test = dd_test))
}

compute_loo_elpd_difference <- function(Ma_fit, Mb_fit) {
  # compute LOO-CV elpd difference
  log_lik_Ma <- extract_log_lik(Ma_fit, merge_chains = FALSE)
  log_lik_Mb <- extract_log_lik(Mb_fit, merge_chains = FALSE)
  r_eff_Ma <- relative_eff(exp(log_lik_Ma))
  r_eff_Mb <- relative_eff(exp(log_lik_Mb))
  loo_Ma <- loo(log_lik_Ma, r_eff = r_eff_Ma, cores = 2)
  loo_Mb <- loo(log_lik_Mb, r_eff = r_eff_Mb, cores = 2)
  loo_elpd_diff <- sum(loo_Ma$pointwise[,"elpd_loo"] 
                       - loo_Mb$pointwise[,"elpd_loo"])
  return(loo_elpd_diff)
}

compute_test_elpd_difference <- function(Ma_fit, Mb_fit) {
  # compute test elpd difference
  test_log_lik_Ma <- extract_log_lik(Ma_fit, 
                                     parameter_name = "log_lik_test")
  test_log_lik_Mb <- extract_log_lik(Mb_fit, 
                                     parameter_name = "log_lik_test")
  elpd_Ma <- elpd(test_log_lik_Ma)
  elpd_Mb <- elpd(test_log_lik_Mb)
  elpd_diff <- sum(elpd_Ma$pointwise[,"elpd"] 
                   - elpd_Mb$pointwise[,"elpd"])
  return(elpd_diff)
}

# compile the models
exec <- cmdstan_model("./stan/K_model_bias.stan")

# define model-fitting method
fit_candidate_model <- function(k, exec, data, n) {
  stan_data <- list(N_train = n,
                    N_test = n,
                    d = 2,
                    x_train = as.matrix(data$train)[, paste0("x", c(0, k))],
                    x_test = as.matrix(data$test)[, paste0("x", c(0, k))],
                    y_train = data$train$y,
                    y_test = data$test$y)
  model_fit <- stanfit(
    exec$sample(data = stan_data,
                chains = 4,
                parallel_chains = 4,
                refresh = 0)
  )
  return(model_fit)
}

one_step <- function(iter, n, K, eps, beta_delta) {
  # simulate data
  data <- simulate_data(n, K, eps, beta_delta)
  
  # fit all models
  baseline_stan_data <- list(N_train = n,
                             N_test = n,
                             d = 1,
                             x_train = as.matrix(data$train$x0),
                             x_test = as.matrix(data$test$x0),
                             y_train = data$train$y,
                             y_test = data$test$y)
  baseline_model <- stanfit(
    exec$sample(data = baseline_stan_data,
                chains = 4,
                parallel_chains = 4,
                refresh = 0)
  )
  fitted_models <- 1:K |> 
    map(\(k) fit_candidate_model(k, exec, data, n))
  
  # compute the difference between all models and the baseline model
  loo_elpd_differences <- fitted_models |>
    map(\(Ma_fit) compute_loo_elpd_difference(Ma_fit, baseline_model))
  test_elpd_differences <- fitted_models |>
    map(\(Ma_fit) compute_test_elpd_difference(Ma_fit, baseline_model))
  
  # build dataframe of results
  out <- cbind(data.frame(sapply(loo_elpd_differences,c)),
               data.frame(sapply(test_elpd_differences,c)))
  names(out) <- c("loo_diff", "test_diff")
  out$model <- 1:K
  out$K <- K
  out$beta <- beta_delta
  return(out)
}

experiment <- function(n, K, num_iters, eps, beta_delta) {
  # concatenate the results
  out <- 1:num_iters |>
    map(\(iter) one_step(iter, n, K, eps, beta_delta)) |>
    bind_rows(.id = "iter")
  return(out)
}

# vary both the number of data observations and the correlation between parameters
n <- 512
eps <- 1
num_iters <- 100
Ks <- c(2, 10, 100)
betas <- seq(0,1,length.out=25)
options <- expand.grid(n=n, K=Ks, num_iters=num_iters, eps=eps, beta=betas)


# run experiment
out <- experiment(n, K, num_iters, eps, beta_delta)