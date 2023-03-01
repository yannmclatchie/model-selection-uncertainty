library(ggplot2)
library(ggdist)
library(patchwork)
library(purrr)
library(dplyr)
library(stringr)
library(posterior)

# load model elppd
( model_df <- readRDS("./data/model_df.rds") )
( p <- ncol(model_df) - 1 )

# set plot params
my_alpha <- 0.05

# compute pointwise elpd difference for the models
( diff_model_df <- model_df - base_fit_elppd )

# compute the elpd different point estimate and standard error for models
( diff_df <- diff_model_df %>%
    reshape2::melt(variable.name = "model", value.name = "diff_elppd") %>%
    group_by(model) %>%
    summarise(
      diff = sum(diff_elppd), 
      diff.se = sd(diff_elppd - diff) * sqrt(n),
      n = n()
    )  )
tail(diff_df)

# plot Gaussian approximation
( gg_diff_dist_analytic <- diff_df %>%
  mutate(
    dist = "norm",
    args = map2(diff, diff.se, list)
  ) %>%
  ggplot(
    aes(
      xdist = dist,
      args = args, 
      colour = model, 
      alpha = model
    )
  ) +
  stat_slab(fill = NA) +
  scale_colour_manual(
    breaks = sprintf("x%1$d", 1:(p + 1)), 
    values = c(rep("black", times = p), c("red"))
  ) +
  scale_alpha_manual(
    breaks = sprintf("x%1$d", 1:(p + 1)), 
    values = c(rep(my_alpha, times = p), 1)
  ) +
  theme_classic() +
  ylab(NULL) +
  xlab(NULL) +
  ggtitle("LOO-CV ELPD difference") +
  theme(axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none") )

# compile hierarchical model
hier_model <- cmdstan_model("./stan/eight_schools_rhs.stan")

# extract data from models
data_list <- list(
  J = nrow(diff_df),
  y = diff_df$diff,
  sigma = diff_df$diff.se,
  hs_df = 3,
  hs_df_global = 1,
  hs_df_slab = 4,
  hs_scale_global = 1,
  hs_scale_slab = 2
)

# fit the hierarchical model
hier_fit <- hier_model$sample(
  data = data_list, 
  seed = 1234, 
  chains = 4, 
  parallel_chains = 4
)

# compute posterior sample statistics
( hier_diff_df <- as_draws_df(hier_fit$draws()) %>%
    select(starts_with("theta")) %>%
    summarise(
      across(
        everything(), 
        list(
          mean = ~mean(.x, na.rm = TRUE),
          sd = ~sd(.x, na.rm = TRUE)
        )
      )
    ) %>%
    reshape2::melt() %>%
    mutate(
      model = str_split(variable, "_", simplify = T)[, 1],
      stat = str_split(variable, "_", simplify = T)[, 2]
    ) %>% 
    select(!variable) %>%
    reshape2::dcast(model ~ stat, value = 'value') %>%
    mutate(model = as.numeric(gsub("\\D", "", model)), dist = "norm") %>%
    arrange(model) %>%
    mutate(model = paste0("x", as.character(model))) )

# plot Gaussian approximation to posterior
( gg_diff_dist_hier <- hier_diff_df %>% ggplot(
    aes(
      xdist = dist_normal(mean, sd),
      colour = model, 
      alpha = model
    )
  ) +
  stat_slab(fill = NA) +
  scale_colour_manual(
    breaks = sprintf("x%1$d", 1:(p + 1)), 
    values = c(rep("black", times = p), c("red"))
  ) +
  scale_alpha_manual(
    breaks = sprintf("x%1$d", 1:(p + 1)), 
    values = c(rep(my_alpha, times = p), 1)
  ) +
  scale_fill_manual(values = c("transparent", "skyblue")) +
  theme_classic() +
  ylab(NULL) +
  xlab(NULL) +
  ggtitle("Eight schools $theta$") +
  theme(axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none") )

# combine plots in patchwork
( gg_elpd_dists <- gg_diff_dist_analytic | gg_diff_dist_hier )

