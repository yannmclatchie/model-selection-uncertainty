library(dplyr)
library(ggplot2)
library(patchwork)
source("./R/aux/aux_plotting.R")

results_df <- read.csv(file = "data/results/jointpred_all.csv")

# make plot

( p <- results_df %>% mutate(is_safe = (loo_elpd_diff >= -4)) %>%
    ggplot(aes(x = loo_elpd_diff, y = elpd_diff, colour = is_safe)) + 
    geom_hline(yintercept=0, colour = "grey") + 
    geom_vline(xintercept=0, colour = "grey") + 
    geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), 
                 data = data.frame(x1 = -140, x2 = 10, y1 = -140, y2 = 10),
                 colour="grey", linetype="dashed") +
    #geom_abline(slope=1, colour="grey", linetype="dashed") +
    geom_vline(xintercept=-4, colour = "red", linetype=2) + 
    geom_point() + 
    #annotate('rect', xmin=0, xmax=Inf, ymin=0, ymax=-Inf, alpha=.3, fill='red') +
    #annotate('rect', xmin=0, xmax=-Inf, ymin=0, ymax=Inf, alpha=.3, fill='red') +
    scale_color_manual(values = c("black", "red")) + 
    ylab("Delta elpd") + 
    xlab("Delta elpd hat") + 
    theme_classic() +
    theme(legend.position="none") )

( p_marginals <- ggExtra::ggMarginal(
  p,
  type = 'density',
  margins = 'y',
  size = 2,
  colour = '#E6000000',
  groupFill = TRUE,
  alpha = 0.7,
  bw = "SJ"
))

# save plot
#save_tikz_plot(p_marginals, width = 5, filename = "./tex/nested-regression-bias.tex")

p_false_positive <- results_df %>% mutate(
  V5 = as.factor(
    ifelse(loo_elpd_diff >= -4 & elpd_diff >= -4, 1, 
              ifelse(loo_elpd_diff >= -4 & elpd_diff < -4, 2,
                     ifelse(loo_elpd_diff < -4 & elpd_diff < -4, 3, 4))) )
  ) %>%
  ggplot(aes(x = loo_elpd_diff, y = elpd_diff, colour = V5)) + 
  geom_hline(yintercept=0, colour = "black") + 
  geom_vline(xintercept=0, colour = "black") +
  geom_vline(xintercept=-4, colour = "grey", linetype=2) + 
  geom_point() + 
  annotate('rect', xmin=-Inf, xmax=-4, ymin=-4, ymax=4, alpha=.2, fill='blue') +
  annotate('rect', xmin=-Inf, xmax=-4, ymin=4, ymax=Inf, alpha=.2, fill='red') +
  scale_x_continuous(limits = c(-25, 5)) + 
  scale_y_continuous(limits = c(-20, 10)) + 
  scale_color_manual(values = c("grey", "grey", "black", "black")) + 
  ylab("Delta elpd") + 
  xlab("Delta elpd hat") + 
  labs(subtitle = "False positives") +
  theme_classic() +
  theme(legend.position="none")
p_false_negative <- results_df %>% mutate(
  V5 = as.factor(
    ifelse(loo_elpd_diff >= -4 & elpd_diff >= -4, 1, 
           ifelse(loo_elpd_diff >= -4 & elpd_diff < -4, 2,
                  ifelse(loo_elpd_diff < -4 & elpd_diff < -4, 3, 4))) )
) %>%
  ggplot(aes(x = loo_elpd_diff, y = elpd_diff, colour = V5)) + 
  geom_hline(yintercept=0, colour = "black") + 
  geom_vline(xintercept=0, colour = "black") +
  geom_vline(xintercept=-4, colour = "grey", linetype=2) + 
  geom_point() + 
  annotate('rect', xmin=-4, xmax=Inf, ymin=-Inf, ymax=-4, alpha=.2, fill='red') +
  annotate('rect', xmin=-4, xmax=Inf, ymin=-4, ymax=4, alpha=.2, fill='blue') +
  scale_x_continuous(limits = c(-25, 5)) + 
  scale_y_continuous(limits = c(-20, 10)) + 
  scale_color_manual(values = c("grey", "grey", "black", "black")) + 
  ylab("") + 
  xlab("Delta elpd hat") + 
  labs(subtitle = "False negatives") +
  theme_classic() +
  theme(legend.position="none")
( p_errors <- p_false_positive + p_false_negative )
#save_tikz_plot(p_errors, width = 5, filename = "./tex/nested-regression-errors.tex")

( p_betas <- results_df %>% mutate(is_safe = (loo_elpd_diff >= -4)) %>%
  ggplot(aes(beta_delta, fill = is_safe)) +
  geom_density(bw = "SJ", colour = NA, alpha = 0.5) +
  annotate(
    "label", label = "Safe",
    x = -2, y = 0.3, size = 3, colour = "black"
  ) +
  annotate(
    "label", label = "Unsafe",
    x = 0.8, y = 1, size = 3, colour = "red"
  ) +
  geom_vline(xintercept=0.4, colour = "grey", linetype=2) + 
  geom_vline(xintercept=-0.4, colour = "grey", linetype=2) + 
  scale_fill_manual(values = c("black", "red")) + 
  ylab("Posterior mean") + 
  xlab("beta Delta") + 
  theme_classic() +
  theme(legend.position="none") )
#save_tikz_plot(p_betas, width = 5, filename = "./tex/nested-regression-betas.tex")

