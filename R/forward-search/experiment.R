library(tidyverse)
library(brms)
library(loo)
library(projpred)

# Config loads dataset configs, hyperparameters and utility functions
source('R/forward-search/config.R')
source('R/forward-search/forward_search.R')
source('R/forward-search/projpred_search.R')

args <- commandArgs(trailingOnly = TRUE)

# read command line arguments
datasets_file <- args[[1]]
n <- args[[2]]
rho <- args[[3]]
prior_name <- args[[4]]
dataset_iter <-  args[[5]]

print(paste0("n: ", n, " rho: ", rho, " prior: ",prior_name))

# ensure correct data types
n <- as.numeric(n)
rho <- as.numeric(rho)
dataset_iter <- as.numeric(dataset_iter)

# read the dataset
datasets <- readRDS(datasets_file)
current_data <- datasets[[dataset_iter]]
df <- current_data$train
df_test <- current_data$test

if (prior_name == "r2d2" || prior_name == "projpred") {
  prior = c(
    prior(R2D2(mean_r2, prec_r2, cons_d2), class = 'b'),
    prior('normal(0, 2.5)', class='Intercept'),
    prior('exponential(1)', class='sigma')
  )
} else if (prior_name == "normal") {
  prior = c(
    prior('normal(0, 1)', class = 'b'),
    prior('normal(0, 2.5)', class='Intercept'),
    prior('exponential(1)', class='sigma')
  )
} else {
  warning("Invalid prior name")
}

reference = brm(y ~ ., 
                data = df, 
                prior = prior, 
                backend = 'cmdstanr',
                #stan_model_args=list(stanc_options = list("O1")),
                #control=list(adapt_delta=0.95),
                silent=2, 
                refresh=0)
ref_elpd_loo = get_elpd_loo(reference)
ref_elpd_test = get_elpd_test(reference, df_test)

base = update(reference, formula=y ~ 1, refresh=0, silent=2)
base$variable = 'base'
base$loo = loo(base)

if (prior_name == "projpred") {
  # Run projpred
  results <- run_projpred_varsel(
    reference, 
    df_test, 
    # Following columns are added as metadata to the final output table
    name=prior_name,
    iter=dataset_iter,
    elpd_loo_ref=ref_elpd_loo,
    elpd_test_ref=ref_elpd_test,
    n=n,
    rho=rho,
    n_test=n_test)
  candidates <- NULL
} else {
  # Run forward selection with LOO-ELPD as the selection criteria
  out = run_forward_selection(
    model = base,
    train_data = df,
    test_data = df_test,
    prior = prior,
    # Following columns are added as metadata to the final output table
    name=prior_name,
    iter=dataset_iter,
    elpd_loo_ref=ref_elpd_loo,
    elpd_test_ref=ref_elpd_test,
    n=n,
    rho=rho,
    n_test=n_test)
  results <- out$results
  candidates <- out$candidates
}

output_file <- paste0("forward_results_",dataset_iter,"_",n,"_",rho,"_",prior_name)
saveRDS(list(forward_search=results,
             candidates=candidates), 
        file = paste0("data/results/forward-search/", output_file, ".RDS"))

print(output_file)
print("done!")
