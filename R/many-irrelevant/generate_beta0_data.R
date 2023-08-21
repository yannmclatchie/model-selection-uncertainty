library(simstudy)
library(bayesflow)

args <- commandArgs(trailingOnly = TRUE)
n <- as.numeric(args[[1]]) # n = 100
n_test <- as.numeric(args[[2]]) # n_test = 1000

simulate_data <- function(rep_id, n, n_test, K, beta_delta = 0) {
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

for (K in seq(2, 106, by=8)) {
    print(paste0("K: ",K))

    # simulate datasets
    datasets <- bayesflow::generate_from_dgp(
        dgp = simulate_data,
        n_datasets = 100,
        n_datasets = 100,
        n = n,
        n_test = n_test,
        K = K,
        beta_delta = 0
    )

    saveRDS(datasets, paste0("data/datasets/all-irrelevant/all_irrelevant_", 
                                "K", K, 
                                "_datasets.RDS"))  
}
