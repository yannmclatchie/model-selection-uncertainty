library(simstudy)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
n <- as.numeric(args[[1]])

simulate_gaussian_data <- function(rep_id, n) {
  # simulate beta from standard normal
  beta_delta <- rnorm(1)
  
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
  dd_train <- genData(n, def)
  dd_test <- genData(n, def)
  
  # return output
  return(list(n = n, train = dd_train, test = dd_test, beta_delta = beta_delta))
}

# simulate datasets
set.seed(1234)
datasets <- bayesflow::generate_from_dgp(
  dgp = simulate_data,
  n_datasets = 1000,
  n = n
)

saveRDS(datasets, paste0("data/datasets/jointpred_n",n,"_datasets.RDS"))