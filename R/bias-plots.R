library(dplyr)
library(tidyr)
library(distributional)
library(ggdist)
library(ggplot2)
library(patchwork)
library(bayesflow)

theme_set(theme_ggdist())

dists_df <- data.frame(
  mean = seq(0, 3, length.out = 10), 
  sd = 1
)

mywidth <- 1

p_dists <-  dists_df %>%
  ggplot(aes(x = mean, ydist = dist_normal(mean, sd))) +
  stat_slab(aes(fill = after_stat(y < 0)), 
            slab_colour = "black",
            slab_linewidth = mywidth) +
  scale_fill_manual(values = c("white", "red")) +
  scale_color_viridis_c() +
  labs(
    y = NULL,
    x = "mean"
  ) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.background = element_blank(),
        legend.position="none")
p_dists

bias_fun <- function(x, mean, sd) return( x * dnorm(x, mean = mean, sd = sd) )
dists_df$bias <- dists_df$mean |> 
  purrr::map_vec(\(mean) integrate(bias_fun, mean = mean, sd = 1, 
                                   lower = -Inf, upper = 0)$value)
p_biases <- dists_df |>
  mutate(bias = abs(bias)) |>
  ggplot(aes(y = bias, x = mean)) +
  geom_point(colour = "red", size = mywidth * 2) +
  geom_line(colour = "red", linewidth = mywidth) +
  labs(
    y = NULL,
    x = "mean"
  ) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.background = element_blank(),
        legend.position="none")
p_biases

( p <- p_dists + p_biases )

# save plot
save_tikz_plot(p, width = 5, height = 3, filename = "./tex/gaussian-bias.tex")
