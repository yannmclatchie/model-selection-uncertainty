library(simstudy)
library(bayesflow)

args <- commandArgs(trailingOnly = TRUE)
n <- as.numeric(args[[1]]) # n = 100
n_test <- as.numeric(args[[2]]) # n_test = 1000

simulate_data <- function(rep_id, n, n_test, K, beta_delta) {
  # compute the SNR to keep marginal variance fixed
  sigma2 <- 1 - beta_delta^2
  
  # define the DGP 
  def <- defData(varname = "x0", formula = "1")
  def <- defRepeat(def, nVars = K - 1, prefix = "x", formula = "0",
                   variance = "1", dist = "normal")
  def <- defData(def, "y", formula = "1 * x0 + ..beta_delta * x1", 
                 variance = "..sigma2", dist = "normal")
  
  # generate the data
  dd_train <- genData(n, def)
  dd_test <- genData(n_test, def)
  
  # return output
  return(list(train = dd_train, 
              test = dd_test, 
              rep_id = rep_id,
              n = n, 
              n_test = n_test,
              K = K, 
              sigma2 = sigma2, 
              snr = beta_delta^2 / sigma2,
              beta_delta = beta_delta))
}

# define values of beta_delta
beta_deltas <- round(10^(seq(-3, -0.001, length.out = 15)), 3)

for (K in c(2, 10, 100)) {
  
  for (beta_delta in beta_deltas) {
    print(paste0("K: ",K," beta: ",beta_delta))
    
    # simulate datasets
    set.seed(1234)
    datasets <- bayesflow::generate_from_dgp(
      dgp = simulate_data,
      n_datasets = 100,
      n = n,
      n_test = n_test,
      K = K,
      beta_delta = beta_delta
    )

    saveRDS(datasets, paste0("data/datasets/many-irrelevant/many_models_", 
                             "K", K, 
                             "_beta", format(beta_delta, nsmall = 2), 
                             "_datasets.RDS"))
  }
  
}
