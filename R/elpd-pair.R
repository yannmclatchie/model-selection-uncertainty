library(dplyr)
library(ggplot2)
library(ggExtra)
library(loo)
library(brms)

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

# compute elpds
elpd_nopool <- loo(posterior_nopool)$pointwise[,"elpd_loo"]
elpd_pool <- loo(posterior_draws)$pointwise[,"elpd_loo"]

llk_pool <- log_lik(posterior_draws)
elpd_pool <- elpd(llk_pool)$pointwise[,"elpd"]
sum(elpd_pool)
sqrt(length(elpd_pool) * var(elpd_pool))
elpd(llk_pool)

# perform standard loo_compare
loo_compare(loo(posterior_nopool), loo(posterior_draws))

## plot the elpds
## --------------

p <- data.frame(elpd1 = elpd_nopool, elpd2 = elpd_pool) %>%
  ggplot() +
  geom_point(aes(x = elpd1, y = elpd2)) +
  annotate(
    "segment", 
    x = min(elpd_nopool, elpd_pool), 
    xend = 0, 
    y = min(elpd_nopool, elpd_pool), 
    yend = 0,
    colour = "blue", 
    linetype = "dashed"
  ) +
  theme_bw()
p <- ggMarginal(p, type = "density")
p
ggsave(
  plot = p, 
  filename = "./img/elpd-pairs.pdf"
)
