library(ggplot2)
library(geomtextpath)
library(bayesflow)

diff.se <- 2

pbm <- function(diff) { 1 / (1 + exp(-diff)) }
pbm_plus <- function(diff, diff.se = 2) { 
  my_fun <- function (z) { dnorm(z, mean = 0, sd = diff.se) / (1 + exp(-diff - z)) }
  out <- integrate(my_fun, -Inf, Inf)
  return ( out$value )
}
elpd_diff_prob <- function(diff) { pnorm(0, mean = diff, sd = diff.se, lower.tail = F) }


t <- seq(0, 10, length.out = 1000)
lookup <- c("pseudo-BMA+" = "pbm_plus", "pseudo-BMA" = "pbm", "Normal approx." = "diff_prob")

( p <- data.frame(
    pbm = pbm(t),
    pbm_plus = t |> map_dbl(\(diff) pbm_plus(diff, diff.se = diff.se)),
    diff_prob = elpd_diff_prob(t)
  ) %>% rename(all_of(lookup)) %>%
  reshape2::melt() %>%
  mutate(t = rep(t, 3)) %>%
  ggplot(aes(x = t, y = value, colour = variable, label = variable)) +
  geom_labelpath(aes(hjust = variable), linewidth = 1.5, size = 3) +
  geom_vline(xintercept = 4, linetype = 2, colour = "grey") +
  scale_x_continuous(name="Delta elpdHatPlain", breaks=c(0, 4)) +
  ylab("Evidence for more complex model") +
  theme_classic() + 
  scale_colour_manual(values = c("red", "darkgrey", "black")) +
  #scale_hjust_discrete(range = c(0.8,0,0.5)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.background = element_blank(),
        legend.position="none") )

save_tikz_plot(p, width = 5, filename = "./tex/elpd-diff-pbma.tex")
