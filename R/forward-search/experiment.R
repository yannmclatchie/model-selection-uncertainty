library(tidyverse)
library(brms)
library(loo)

args <- commandArgs(trailingOnly = TRUE)

datasets_file <- args[[1]]
dataset_iter <- args[[2]]
dataset_iter <- as.numeric(dataset_iter)

datasets <- readRDS(datasets_file)

current_data <- datasets[[dataset_iter]]
df <- current_data$train
df_test <- current_data$test
 
# Config loads dataset configs, hyperparameters and utility functions
source('R/forward-search/config.R')
source('R/forward-search/forward_search.R')

prior = c(
  prior(R2D2(mean_r2, prec_r2, cons_d2), class = 'b'),
  prior('normal(0, 2.5)', class='Intercept'),
  prior('exponential(1)', class='sigma')
)

prior_normal = c(
  prior('normal(0, 1)', class = 'b'),
  prior('normal(0, 2.5)', class='Intercept'),
  prior('exponential(1)', class='sigma')
)

reference = brm(y ~ ., 
                data = df, 
                prior = prior, 
                backend = 'cmdstanr',
                stan_model_args=list(stanc_options = list("O1")),
                control=list(adapt_delta=0.95),
                silent=2, 
                refresh=0)
ref_elpd_loo = get_elpd_loo(reference)
ref_elpd_test = get_elpd_test(reference, df_test)

base = update(reference, formula=y ~ 1, refresh=0, silent=2)
base$variable = 'base'
base$loo = loo(base)

# Run forward selection with LOO-ELPD as the selection criteria
out = run_forward_selection(
  model = base,
  train_data = df,
  test_data = df_test,
  prior = prior,
  # Following columns are added as metadata to the final output table
  name='R2D2 prior',
  iter=paste(job_id, array_id, sep='_'),
  elpd_loo_ref=ref_elpd_loo,
  elpd_test_ref=ref_elpd_test,
  n=n,
  rho=rho,
  n_test=n_test)

# Run forward selection with LOO-ELPD as the selection criteria
out_normal = run_forward_selection(
  model = base,
  train_data = df,
  test_data = df_test,
  prior = prior_normal,
  # Following columns are added as metadata to the final output table
  name='Normal prior',
  iter=paste(job_id, array_id, sep='_'),
  elpd_loo_ref=ref_elpd_loo,
  elpd_test_ref=ref_elpd_test,
  n=n,
  rho=rho,
  n_test=n_test)


output_file <- paste0("forward_search_results_iter", dataset_iter)
saveRDS(list(forward_search=out$results,
             candidates=out$candidates,
             forward_search_normal=out_normal$results), 
        file = paste0("data/results/", output_file, ".RDS"))
