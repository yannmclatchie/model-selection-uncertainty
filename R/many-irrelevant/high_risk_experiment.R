library(cmdstanr)
library(purrr)
library(loo)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

datasets_file <- args[[1]]
dataset_iter <- args[[2]]
dataset_iter <- as.numeric(dataset_iter)
exec <- args[[3]]

datasets <- readRDS(datasets_file)

current_data <- datasets[[dataset_iter]]

compute_loo_elpd_difference <- function(Ma_fit, Mb_fit, n) {
  # compute LOO-CV elpd difference
  log_lik_Ma <- Ma_fit$draws("log_lik")
  log_lik_Mb <- Mb_fit$draws("log_lik")
  r_eff_Ma <- relative_eff(exp(log_lik_Ma))
  r_eff_Mb <- relative_eff(exp(log_lik_Mb))
  loo_Ma <- loo(log_lik_Ma, r_eff = r_eff_Ma, cores = 2)
  loo_Mb <- loo(log_lik_Mb, r_eff = r_eff_Mb, cores = 2)
  loo_elpd_diff <- sum(loo_Ma$pointwise[,"elpd_loo"] 
                       - loo_Mb$pointwise[,"elpd_loo"])
  loo_elpd_diff_se  <- sd(loo_Ma$pointwise[, 'elpd_loo'] - 
                          loo_Mb$pointwise[, 'elpd_loo']) * sqrt(n)
  return(list(loo_elpd_diff = loo_elpd_diff, loo_elpd_diff_sd = loo_elpd_diff_se))
}

compute_test_elpd_difference <- function(Ma_fit, Mb_fit) {
  # compute test elpd difference
  test_log_lik_Ma <- Ma_fit$draws("log_lik_test")
  test_log_lik_Mb <- Mb_fit$draws("log_lik_test")
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
  model_fit <- exec$sample(data = stan_data,
                           chains = 4,
                           parallel_chains = 4,
                           refresh = 0)
  return(model_fit)
}

fit_all_models <- function(exec_file, data) {
  n <- as.numeric(data$n)
  K <- as.numeric(data$K)
  iter <- data$rep_id
  
  # build model
  exec <- cmdstan_model(exe_file = exec_file)
  
  # fit all models
  print("fitting baseline model ...")
  baseline_stan_data <- list(N_train = n,
                             N_test = n,
                             d = 1,
                             x_train = as.matrix(data$train$x0),
                             x_test = as.matrix(data$test$x0),
                             y_train = data$train$y,
                             y_test = data$test$y)
  baseline_model <- exec$sample(data = baseline_stan_data,
                                chains = 4,
                                parallel_chains = 4,
                                refresh = 0)
  print("fitting candidate models ...")
  fitted_models <- 1:K |> 
    map(\(k) fit_candidate_model(k, exec, data, n))
  print("done.")
  
  # compute the difference between all models and the baseline model
  out <- fitted_models |>
    map(\(Ma_fit) compute_loo_elpd_difference(Ma_fit, baseline_model, n)) |>
    bind_rows()

  # build dataframe of results
  out$K <- K
  out$iter <- iter
  return(out)
}

# fit all models and compute stats
out <- fit_all_models(exec, current_data)

# save results
K <- as.numeric(current_data$K)
beta_delta <- as.numeric(current_data$beta_delta)
output_file <- paste0("all_irrelevant_results_","K",K,
                      "_iter", dataset_iter)
write.csv(out, file = paste0("data/results/all-irrelevant/", output_file, ".csv"), 
          row.names = FALSE)

