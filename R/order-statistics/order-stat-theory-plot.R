library(ggplot2)
library(dplyr)
library(bayesflow)

# define order statistic heuristic data
order_stat_heuristic <- function(k, c = 2.8) {
  qnorm(p = 1 - 1 / (k * 2), mean = 0, sd = c)
}

# define the observed data
x <- 2:20

# define plotting data
order_stat_df <- data.frame(K = x, 
                            # we multiply the order stat by 1.5 to estimate selection bias
                            y = x |> purrr::map_dbl(\(k) order_stat_heuristic(k, c = 1) * 1.5))
(p <- ggplot(data = order_stat_df, aes(x = K, y = y)) + 
    geom_line(colour = "red") +
    geom_point(colour = "red") +
    ylab("Expected\nselection-induced bias") + 
    xlab("$K$") + 
    scale_x_continuous(trans='log2') +
    #scale_y_continuous(breaks=c(0,5,10,15,20),limits=c(4,12)) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank()) )

# save plot
save_tikz_plot(p, width = 5, height = 2, filename = "./tex/order-stat-theory.tex")
