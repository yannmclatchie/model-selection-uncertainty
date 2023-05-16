library(cmdstanr)
library(argparse)

parser <- ArgumentParser(description = "compile models")
parser$add_argument("stan_path", metavar = "O", type = "character", help = "stan output file path")

args <- parser$parse_args()

stan_path <- args$stan_path

models <- list(
  K_model_bias = cmdstan_model(
    "/scratch/work/mclatcy1/model-selection-uncertainty/stan/K_model_bias.stan", 
    # dir = stan_path, 
    # stanc_options = list("O1"),
    force_recompile = TRUE),
  joint_predictive_lin_reg = cmdstan_model(
    "/scratch/work/mclatcy1/model-selection-uncertainty/stan/joint_predictive_lin_reg.stan", 
    # dir = stan_path, 
    # stanc_options = list("O1"),
    force_recompile = TRUE),
  hier_elpd = cmdstan_model(
    "/scratch/work/mclatcy1/model-selection-uncertainty/stan/hier-elpd.stan", 
    # dir = stan_path, 
    # stanc_options = list("O1"),
    force_recompile = TRUE),
  eight_schools = cmdstan_model(
    "/scratch/work/mclatcy1/model-selection-uncertainty/stan/eight_schools.stan", 
    # dir = stan_path, 
    # stanc_options = list("O1"),
    force_recompile = TRUE),
  eight_schools_rhs = cmdstan_model(
    "/scratch/work/mclatcy1/model-selection-uncertainty/stan/eight_schools_rhs.stan", 
    # dir = stan_path, 
    # stanc_options = list("O1"),
    force_recompile = TRUE)
)

print("done!")
