library(simstudy)
library(bayesflow)

args <- commandArgs(trailingOnly = TRUE)
n <- as.numeric(args[[1]])
eps <- as.numeric(args[[2]])

simulate_data <- function(rep_id, n, K, eps, beta_delta = 0) {
  # define the DGP 
  def <- defData(varname = "x0", formula = "10", variance = "10", dist = "normal")
  def <- defRepeat(def, nVars = K, prefix = "x", formula = "..eps",
                   variance = "10", dist = "normal")
  def <- defData(def, "yInlier", formula = "1 + x0 + ..beta_delta * x1", 
                 variance = "..eps", dist = "normal")
  def <- defData(def, varname = "yOutlier", formula = "50", 
                 variance = "..eps", dist = "normal")
  def <- defData(def, varname = "y", 
                 formula = "yInlier | .8 + yOutlier | .2", 
                 dist = "mixture")

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

for (K in seq(2, 106, by=8)) {
    print(paste0("K: ",K))

    # simulate datasets
    set.seed(1234)
    datasets <- bayesflow::generate_from_dgp(
        dgp = simulate_data,
        n_datasets = 100,
        n = n,
        K = K,
        eps = eps,
        beta_delta = 0
    )

    saveRDS(datasets, paste0("data/datasets/all-irrelevant/all_irrelevant_", 
                                "K", K, 
                                "_datasets.RDS"))  
}