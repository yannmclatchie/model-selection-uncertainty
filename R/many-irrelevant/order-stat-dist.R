library(ggplot2)

# define regime
p <- 101
mu <- 0
sd <- 2.3

# with ppoints 
max_ppoint <- max(ppoints(p, a = 0.5))
( qnorm(p = max_ppoint, mean = mu, sd = sd) )

# manually
( qnorm(p = 1 - 1 / (p * 2), mean = mu, sd = sd) )

# order statistics from uniform
p <- 101
n <- p
beta_samples <- rbeta(n = n, shape1 = p, shape2 = n + 1 - p)
# make density plot
( gg_order_stat_dist <- ggplot(
  data = data.frame(samples = qnorm(p = beta_samples, mean = mu, sd = sd)),
  aes(x = samples)
  ) +
  stat_density(fill = NA, colour = "black", outline.type = "upper") +
  geom_vline(
    xintercept = qnorm(p = p / (n + 1), mean = mu, sd = sd),
    colour = "red",
    linetype = "dashed"
  ) +
  annotate(
    "text", 
    x = round(qnorm(p = p / (n + 1), mean = mu, sd = sd), 2) * 1.1, 
    y = 0.05, 
    label = "Theoretic\nmean",
    colour = "red"
  ) +
  theme_classic() +
  ylab(NULL) +
  xlab("Maximum order statistic") +
  theme(axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) )
  
