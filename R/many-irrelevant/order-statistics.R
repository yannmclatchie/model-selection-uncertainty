library(ggplot2)
library(ggdist)
library(purrr)
library(dplyr)

# load model elppd
( model_df <- readRDS("./data/model_df.rds") )

# compute pointwise elpd difference for the models
( diff_model_df <- model_df - base_fit_elppd )

# plot elpd densities
( gg_diff_dist <- model_df %>%
    reshape2::melt(variable.name = "model", value.name = "elppd") %>%
    ggplot(aes(
      x = elppd, 
      model = model, 
      colour = model, 
      group = model, 
      alpha = model
    )
    ) +
    stat_density(
      geom = "line",
      position = "identity",
      bw = "nrd"
    ) +
    scale_colour_manual(
      breaks = sprintf("x%1$d", 1:(p + 1)), 
      values = c(rep("black", times = p), c("red"))
    ) +
    scale_alpha_manual(
      breaks = sprintf("x%1$d", 1:(p + 1)), 
      values = c(rep(0.2, times = p), 1)
    ) +
    theme_classic() +
    ylab(NULL) +
    xlab("ELPD") +
    theme(axis.line.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "none") )

# plot elpd difference densities
( gg_diff_elpd_dist <- diff_model_df %>%
  reshape2::melt(variable.name = "model", value.name = "diff_elppd") %>%
  ggplot(aes(
    x = diff_elppd, 
    model = model, 
    colour = model, 
    group = model, 
    alpha = model
    )
  ) +
    stat_density(
      geom = "line",
      position = "identity",
      bw = "bcv"
    ) +
    scale_colour_manual(
      breaks = sprintf("x%1$d", 1:(p + 1)), 
      values = c(rep("black", times = p), c("red"))
    ) +
    scale_alpha_manual(
      breaks = sprintf("x%1$d", 1:(p + 1)), 
      values = c(rep(0.2, times = p), 1)
    ) +
    theme_classic() +
    ylab(NULL) +
    xlab("ELPD difference") +
    theme(axis.line.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "none") )

# compute the elpd different point estimate and standard error for models
( diff_df <- diff_model_df %>%
  reshape2::melt(variable.name = "model", value.name = "diff_elppd") %>%
  group_by(model) %>%
  summarise(
    diff = sum(diff_elppd), 
    diff.se = sd(diff_elppd - diff) * sqrt(n),
    n = n()
  )  )

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
