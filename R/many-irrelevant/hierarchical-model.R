library(ggplot2)
library(ggdist)
library(purrr)
library(dplyr)
library(stringr)
library(posterior)

# load model elppd
( model_df <- readRDS("./data/model_df.rds") )
( p <- ncol(model_df) - 1 )

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
diff_df %>%
  mutate(
    dist = "norm",
    args = map2(diff, diff.se, list)
  ) %>%
  ggplot(
    aes(
      y = dist, 
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
    values = c(rep(0.05, times = p), 1)
  ) +
  theme_classic() +
  ylab(NULL) +
  xlab("Model ELPD difference\n(Gaussian approximation)") +
  theme(axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none")

# compile hierarchical model
hier_model <- cmdstan_model("./stan/eight_schools_rhs.stan")

# extract data from models
data_list <- list(
  J = nrow(diff_df),
  y = diff_df$diff,
  sigma = diff_df$diff.se,
  hs_df = 3,
  hs_df_global = 3,
  hs_df_slab = 3,
  hs_scale_global = 10,
  hs_scale_slab = 10
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
hier_diff_df %>% ggplot(
  aes(
    y = dist, 
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
    values = c(rep(0.05, times = p), 1)
  ) +
  theme_classic() +
  ylab(NULL) +
  xlab("Eight schools theta\n(Gaussian approximation)") +
  theme(axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none")
