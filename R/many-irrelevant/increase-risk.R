library(simstudy)
library(cmdstanr)
library(brms)
library(purrr)
library(loo)
stanfit <- function(fit) rstan::read_stan_csv(fit$output_files())


simulate_data <- function(n, K, eps) {
  # define the DGP 
  def <- defRepeat(nVars = 3, prefix = "w", formula = "0",
                   variance = "1", dist = "normal")
  def <- defRepeat(def, nVars = K, prefix = "x", formula = "0",
                   variance = "1", dist = "normal")
  def <- defData(def, "y", formula = "1 * w1 + 0.5 * w2 - 0.3 * w3", 
                 variance = "..eps", dist = "normal")
  
  # generate the data
  dd_train <- genData(n, def)
  dd_test <- genData(n, def)
  
  # return output
  return(list(train = dd_train, 
              test = dd_test, 
              n = n, 
              K = K, 
              eps = eps))
}
data <- simulate_data(n = 10, K = 100, eps = 1)

n <- as.numeric(data$n)
K <- as.numeric(data$K)

exec_file <- "stan/K_model_bias"
exec <- cmdstan_model(exe_file = exec_file)

baseline_stan_data <- list(N_train = n,
                           N_test = n,
                           d = 2,
                           x_train = as.matrix(data$train)[, paste0("x", c(1,2))],
                           x_test = as.matrix(data$train)[, paste0("x", c(1,2))],
                           y_train = data$train$y,
                           y_test = data$test$y)
baseline_model <- stanfit(
  exec$sample(data = baseline_stan_data,
              chains = 4,
              parallel_chains = 4,
              refresh = 0)
)

d <- 10
stan_data <- list(N_train = n,
                  N_test = n,
                  d = d,
                  x_train = as.matrix(data$train)[, paste0("x", 1:d)],
                  x_test = as.matrix(data$test)[, paste0("x", 1:d)],
                  y_train = data$train$y,
                  y_test = data$test$y)
model_fit <- stanfit(
  exec$sample(data = stan_data,
              chains = 4,
              parallel_chains = 4,
              refresh = 0)
)

log_lik_Ma <- extract_log_lik(model_fit, merge_chains = FALSE)
log_lik_Mb <- extract_log_lik(baseline_model, merge_chains = FALSE)
r_eff_Ma <- relative_eff(exp(log_lik_Ma))
r_eff_Mb <- relative_eff(exp(log_lik_Mb))
loo_Ma <- loo(log_lik_Ma, r_eff = r_eff_Ma, cores = 2)
loo_Mb <- loo(log_lik_Mb, r_eff = r_eff_Mb, cores = 2)
loo_compare(loo_Ma, loo_Mb)




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
  baseline_model <- stanfit(
    exec$sample(data = baseline_stan_data,
                chains = 4,
                parallel_chains = 4,
                refresh = 0)
  )
  print("fitting candidate models ...")
  fitted_models <- 1:K |> 
    map(\(k) fit_candidate_model(k, exec, data, n))
  print("done.")
  
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
  return(out)
}

# define model exec
exec <- "stan/increase_risk"

# fit all models and compute stats
out <- fit_all_models(exec, current_data)