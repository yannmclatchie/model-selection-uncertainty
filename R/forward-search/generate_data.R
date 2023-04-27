library(bayesflow)

# Import data generating settings
source('R/forward-search/config.R')

# Generate data according to Piironen and Vehtari 2017
generate_data = function(rep_id, n, rho) {
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

datasets <- bayesflow::generate_from_dgp(
    dgp = generate_data,
    n_datasets = 50,
    n = n,
    rho = rho
)

saveRDS(datasets, paste0("data/datasets/forward_search_dataset.RDS"))
