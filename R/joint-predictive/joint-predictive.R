library(cmdstanr)
library(simstudy)
library(loo)
library(ggplot2)
library(dplyr)
stanfit <- function(fit) rstan::read_stan_csv(fit$output_files())

# define linear regression models
model <- "data {
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
}"

# compile the models
stan_file <- write_stan_file(model)
exec <- cmdstan_model(stan_file)

simulate_gaussian_data <- function(n_train, n_test, beta_delta) {
  # define the DGP 
  def <- defData(varname = "x1", formula = "1")
  def <- defData(def, varname = "x2", formula = "0", 
                 variance = "1", dist = "normal")
  def <- defData(def, varname = "x3", formula = "0", 
                 variance = "1", dist = "normal")
  ( def <- defData(def, "y", 
                   formula = paste0("0 * x1 + 1 * x2 + ..beta_delta * x3"), 
                   variance = "1", dist = "normal") ) 
  
  # generate the data
  dd_train <- genData(n_train, def)
  dd_test <- genData(n_test, def)
  
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

experiment <- function(exec, n_test, n_train, beta_delta) {
  # simulate the data
  data <- simulate_gaussian_data(n_train, n_test, beta_delta)
  
  # define stan data
  stan_data_Ma <- list(
    N_train = n_train,
    N_test = n_test,
    d = 2,
    x_train = data$train[,c("x1","x2")],
    y_train = data$train$y,
    x_test = data$test[,c("x1","x2")],
    y_test = data$test$y
  )
  stan_data_Mb <- list(
    N_train = n_train,
    N_test = n_test,
    d = 3,
    x_train = data$train[,c("x1","x2","x3")],
    y_train = data$train$y,
    x_test = data$test[,c("x1","x2","x3")],
    y_test = data$test$y
  )
  
  # compile and fit the models
  Ma_fit <- stanfit(
    exec$sample(data = stan_data_Ma,
                chains = 4,
                parallel_chains = 4,
                refresh = 0)
  )
  Mb_fit <- stanfit(
    exec$sample(data = stan_data_Mb,
                chains = 4,
                parallel_chains = 4,
                refresh = 0)
  )
  
  # compute elpd difference
  loo_elpd_diff <- compute_loo_elpd_difference(Ma_fit, Mb_fit)
  test_elpd_diff <- compute_test_elpd_difference(Ma_fit, Mb_fit)
  
  # save results
  out <- data.frame(beta_delta = beta_delta, 
                    loo_elpd_diff = loo_elpd_diff, 
                    elpd_diff = test_elpd_diff)
  return(out)
}

# define experiment parameters
n_train <- 512
n_test <- n_train
num_iters <- 100
betas <- c(0, 0.5, 1)

# initialise results matrix
results_df <- data.frame(
  matrix(
    vector(), 
    0, 
    3, 
    dimnames=list(c(), c("beta_delta", "loo_elpd_diff", "elpd_diff"))
    )
  )

# run experiment
for (beta_delta in betas) {
  for (iter in 1:num_iters) {
    experiment_df <- experiment(exec, n_test, n_train, beta_delta)
    results_df <- rbind(results_df, experiment_df)
  }
}

# save the data
saveRDS(results_df, file = "data/nested_bias.rds")
results_df <- readRDS(file = "data/nested_bias.rds")

# make plot
cols <- c("0" = "black", "0.5" = "blue", "1" = "red")
( p <- results_df %>% mutate_at("beta_delta", funs(as.factor)) %>%
  ggplot(aes(x = loo_elpd_diff, y = elpd_diff, colour = beta_delta)) + 
  geom_hline(yintercept=0, colour = "grey") + 
  geom_vline(xintercept=0, colour = "grey") + 
  geom_abline(slope=1, colour="grey", linetype="dashed") +
  annotate('rect', xmin=0, xmax=Inf, ymin=0, ymax=-Inf, alpha=.2, fill='red') +
  annotate('rect', xmin=0, xmax=-Inf, ymin=0, ymax=Inf, alpha=.2, fill='red') +
  geom_point() + 
  annotate(
    "label", label = "beta delta = 0",
    x = -10, y = 5, size = 3, colour = "black"
  ) +
  annotate(
    "label", label = "beta delta = 0.5",
    x = -30, y = -10, size = 3, colour = "blue"
  ) +
  annotate(
    "label", label = "beta delta = 1",
    x = -30, y = -60, size = 3, colour = "red"
  ) +
  scale_color_manual(values = cols) + 
  ylab("Delta elpd") + 
  xlab("Delta elpd hat") + 
  coord_cartesian(xlim =c(-70, 5), ylim = c(-70, 5)) +
  theme_classic() +
  theme(legend.position="none") )

# save plot
source("./R/aux/aux_plotting.R")
save_tikz_plot(p, width = 5, filename = "./tex/nested-regression-bias.tex")

