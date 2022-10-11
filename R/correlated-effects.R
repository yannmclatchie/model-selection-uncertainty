library(dplyr)
library(ggplot2)
library(ggExtra)
library(loo)
library(brms)
library(rstan)
library(bayesplot)
library(priorsense)

set.seed(300416)
setwd("~/Desktop/")

## data processing
## ---------------

acti_data <- read.csv("./projpred-workflow/data/monkeys/activity_data.csv") 
activity_2mins <- acti_data %>%
  filter(obs<9) %>% group_by(subj_id, Day) %>%
  summarize(total=sum(Activity), 
            active_bins = sum(Activity > 0), 
            age = min(age)) %>%
  rename(monkey = subj_id, day = Day) %>%
  ungroup()

age_centre <- mean(activity_2mins$age)
age_scale <- diff(range(activity_2mins$age))/2
active_bins_centre <- 4

activity_2mins_scaled <- activity_2mins %>%
  mutate(monkey = factor(monkey),
         day = factor(day),
         age_centred = (age - age_centre)/age_scale,
         active_bins_scaled = (active_bins - active_bins_centre)/4)

## model fitting
## -------------

# define the priors
priors_lm <-  prior(normal(0,1), class = "b") +
  prior(normal(0, 0.2), coef = "age_centred") + 
  prior(normal(0,0.2), coef = "age_centred:day2") +
  prior(normal(0, 1), coef = "day2") +
  prior(normal(0,1), class = "Intercept") +
  prior(normal(0,1), class = "sigma")
priors <- prior(normal(0, 0.2), coef = "age_centred") + 
  prior(normal(0,0.2), coef = "age_centred:day2") +
  prior(normal(0, 1), coef = "day2") +
  prior(normal(0,1), class = "sigma") +
  prior(exponential(1), class = sd) + # tau
  prior(normal(0,1), class = "Intercept")

# fit the models
posterior_nopool <- brm(
  active_bins_scaled ~ age_centred * day + monkey, 
  data = activity_2mins_scaled,
  prior = priors_lm)
posterior_draws <- brm(
  active_bins_scaled ~ age_centred*day + (1 | monkey), 
  data = activity_2mins_scaled,
  prior = priors
)

# extract elpds
extract_elpd <- function(fit) {
  llk <- log_lik(fit)
  elpd <- elpd(llk)$pointwise[,"elpd"]
  return (elpd)
}
elpd_pool <- extract_elpd(posterior_draws)
elpd_nopool <- extract_elpd(posterior_nopool)
elpds <- cbind(elpd_pool, elpd_nopool)

# fit the meta-model
N <- length(activity_2mins_scaled$active_bins_scaled)
data <- list(
  N = N,
  K = 2,
  y = elpds
)
model <- stan_model(file = "correlated-beta.stan")
fit <- sampling(model, data = data, seed = 300416)

# visualise MCMC diagnostics
rhats <- rhat(fit)
effect_params <- c("common_elpd", "beta[1]", "beta[2]")
mcmc_areas(fit, pars = effect_params)
mcmc_rhat(rhats)

# compute the LOO-CV ELPD of the meta-model
log_lik <- extract_log_lik(fit, merge_chains = FALSE)
r_eff <- relative_eff(exp(log_lik), cores = 2) 
loo <- loo(log_lik, r_eff = r_eff, cores = 2)
print(loo)

# prior sensitivity analysis of the meta-model
pss <- powerscale_sequence(fit)
powerscale_plot_ecdf(pss, variables = effect_params) +
  theme_bw() + 
  theme(
    title = element_blank(), 
    panel.grid.minor = element_blank()
  )
