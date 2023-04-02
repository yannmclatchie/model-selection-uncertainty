library(simstudy)
library(bayesflow)

args <- commandArgs(trailingOnly = TRUE)
n <- as.numeric(args[[1]])
eps <- as.numeric(args[[2]])

simulate_data <- function(rep_id, n, K, eps, beta_delta) {
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
  return(list(train = dd_train, 
              test = dd_test, 
              rep_id = rep_id,
              n = n, 
              K = K, 
              eps = eps, 
              beta_delta = beta_delta))
}

for (K in c(2, 10, 100)) {
  
  #for (beta_delta in seq(0, 1, by = 0.05)) {
  for (beta_delta in c(0.0, 1.0)) {  
    # simulate datasets
    set.seed(1234)
    datasets <- bayesflow::generate_from_dgp(
      dgp = simulate_data,
      n_datasets = 5,
      n = n,
      K = K,
      eps = eps,
      beta_delta = beta_delta
    )

    saveRDS(datasets, paste0("data/datasets/many_models_", 
                             "K", K, 
                             "_beta", format(beta_delta, nsmall = 1), 
                             "_datasets.RDS"))
  }
  
}
