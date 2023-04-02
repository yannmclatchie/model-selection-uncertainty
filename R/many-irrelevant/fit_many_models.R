library(cmdstanr)
stanfit <- function(fit) rstan::read_stan_csv(fit$output_files())

args <- commandArgs(trailingOnly = TRUE)

datsets_file <- args[[1]]
dataset_iter <- args[[2]]
exec <- args[[3]]

datasets <- readRDS(datasets_file)

current_data <- datasets[[dataset_iter]]


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

fit_all_models <- function(exec_file, data) {
  n <- as.numeric(data$n)
  K <- as.numeric(data$K)
  beta_delta <- data$beta_delta
  iter <- data$rep_id
  
  # build model
  exec <- cmdstan_model(exe_file = exec_file)

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
  names(out) <- c("loo_elpd_diff", "test_elpd_diff")
  out$model <- 1:K
  out$K <- K
  out$iter <- iter
  out$beta <- beta_delta
  return(out)
}

# fit all models and compute stats
out <- fit_all_models(exec, current_data)

# save results
K <- as.numeric(current_data$K)
beta_delta <- as.numeric(current_data$beta_delta)
output_file <- paste0("many_models_results_","K",K,"_beta", beta_delta, 
                      "_iter", dataset_iter)
saveRDS(out, file = paste0("data/results/", output_file, ".RDS"))


n <- 512
eps <- 1
num_iters <- 1
K <- 2
beta_delta <- 1
curr_data <- simulate_data(n, K, eps, beta_delta)
curr_data$n <- as.character(n)
curr_data$K <- as.character(K)
exec_file <-  "./stan/K_model_bias"
fit_all_models(exec_file, curr_data)
