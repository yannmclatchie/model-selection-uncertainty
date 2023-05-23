library(bayesflow)

args <- commandArgs(trailingOnly = TRUE)

num_iters <- as.numeric(args[[1]])

# Import data generating settings
source('R/forward-search/config.R')

# Generate data according to Piironen and Vehtari 2017
generate_data = function(n, rho) {
  # number of parameters in the model
  p <- 100
  
  # define xi given rho
  if (rho == 0){
    xi <- 0.59
  } else if (rho == 0.5){
    xi <- 0.34
  } else if (rho == 0.9){
    xi <- 0.28
  } else {
    warning("Invalid rho provided.")
  }
  
  # build correlation matrix
  block_size <- 5
  num_matrices <- p / block_size
  listOfMatrices <- vector("list", num_matrices)
  for (i in 1:num_matrices) {
    listOfMatrices[[i]] <- matrix(
      rep(rho, block_size * block_size), nrow=block_size
    )
  }
  R <- matrix(Matrix::bdiag(listOfMatrices), nrow=p)
  diag(R) <- 1
  
  # build associated covariate weights
  w1 <- rep(xi, block_size)
  w2 <- rep(xi * 0.5, block_size)
  w3 <- rep(xi * 0.25, block_size)
  w4 <- rep(0, (num_matrices - 3) * block_size)
  w <- c(w1, w2, w3, w4)
  
  # define zero mean vector
  sigma <- 1
  mu <- rep(0, p)
  sds <- rep(sigma, p)
  
  # sample data
  X <- MASS::mvrnorm(n=n, mu=mu, Sigma=R)
  y <- X %*% w + sigma * rnorm(n=n)
  data <- data.frame(y = y, X = X)
  
  return(data)
}

generate_train_and_test <- function(rep_id, seed, n, n_test, rho) {
  # Keep test data fixed over different jobs...
  df_test = withr::with_seed(seed, generate_data(n_test, rho))
  # ...but randomize over array jobs
  df = generate_data(n, rho)
  return(list(train = df, test = df_test))
}

for (n in c(100, 200, 400)) {
  for (rho in c(0.0, 0.5, 0.9)) {
    datasets <- bayesflow::generate_from_dgp(
      dgp = generate_train_and_test,
      n_datasets = num_iters,
      seed = seed,
      n = n,
      n_test = n_test,
      rho = rho
    )
    dataset_file <- paste0("data/datasets/forward/forward_dataset_",n,"_",rho,".RDS")
    print(dataset_file)
    saveRDS(datasets, dataset_file)
  }
}
