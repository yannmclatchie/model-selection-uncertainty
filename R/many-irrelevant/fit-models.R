library(cmdstanr)
library(brms)
library(loo)
library(purrr)
library(simstudy)

# define experiment settings
n <- 100
p <- 100

# define the DGP 
def <- defRepeat(nVars = p, prefix = "x", formula = "0",
                 variance = "1", dist = "normal")
def <- defData(def, varname = paste0("x", p+1), formula = "0", 
               variance = "1", dist = "normal")
( def <- defData(def, "y", formula = paste0("x", p+1), 
                 variance = "5", dist = "normal") ) # add lots of noise

# generate the data
dd <- genData(n, def)

# define prior for all models
bprior <- c(prior_string("normal(0,1)", class = "b"))

# fit the base model
base_fit <- brm(
  "y ~ 1", 
  data = dd, 
  chains = 4, 
  cores = 4,
  backend = "cmdstanr", 
  threads = threading(4)
)
base_fit_loo <- loo(base_fit)
base_fit_elppd <- base_fit_loo$pointwise[,"elpd_loo"]

# now fit all misspecified models
fit_model <- function(x) {
  new_formula <- paste0("y ~ x", x)
  fit <- update(
    base_fit, 
    formula = new_formula, 
    newdata = dd,
    prior = bprior,
    chains = 4, 
    cores = 4, 
    backend = "cmdstanr", 
    threads = threading(4))
  fit_loo <- loo(fit)
  fit_elppd <- fit_loo$pointwise[,"elpd_loo"]
  return(fit_elppd)
}
xs <- 1:(p + 1)
names(xs) <- paste0("x", xs)
( model_df <- map_dfr(xs, fit_model) )
saveRDS(model_df, "./data/model_df.rds")
