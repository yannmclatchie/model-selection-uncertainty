library(ggplot2)
library(geomtextpath)
library(dplyr)
library(bayesflow)

## Many irrelevant
## ---------------

# read in results
results <- read.csv(file = "data/results/many_models_all.csv")

# plot results
plotting_df <- results %>% group_by(K, beta, iter) %>%
  mutate(best_loo = max(loo_elpd_diff)) %>%
  group_by(K, beta) %>%
  mutate(mean_best_loo = median(best_loo)) %>%
  right_join(results %>% 
               filter(model==1) %>%
               group_by(K, beta, iter) %>%
               mutate(true_diff = median(test_elpd_diff)) %>%
               group_by(K, beta) %>%
               mutate(mean_true_diff = median(true_diff))) %>%
  select(c("beta", "K", "mean_best_loo", "mean_true_diff")) %>%
  unique() %>%
  reshape2::melt(id = c("beta", "K")) %>%
  mutate_at(c("variable"), funs(recode(., `mean_best_loo`="loo", 
                                       `mean_true_diff`="true")))
ann_text_1 <- data.frame(beta = 0.01,value = 5.5, lab = "Text",variable="loo",
                       K = factor(100,levels = c("2","10","100")))
ann_text_2 <- data.frame(beta = 0.01,value = 0.75,lab = "Text",variable="true",
                         K = factor(100,levels = c("2","10","100")))
scaleFUN <- function(x) sprintf("%.2f", x)
( p <- plotting_df %>% 
    ggplot(aes(x = beta, y = value, colour = variable)) + 
    geom_point() +
    geom_smooth(method = "gam", 
                formula = y ~ s(x, k = 6, bs = "cs", m = 5), 
                se = FALSE) +
    facet_wrap(vars(K), labeller = labeller(K = 
                                              c("2" = "$K = 2$",
                                                "10" = "$K = 10$",
                                                "100" = "$K = 100$")
    )) +
    geom_label(data = ann_text_1,label = "Best LOO-CV",size=3) +
    geom_label(data = ann_text_2,label = "True model",size=3,colour="grey") +
    scale_y_continuous(trans='pseudo_log') +
    scale_x_continuous(trans='log10', labels=scaleFUN) +
    ylab("elpd difference") + 
    xlab("$beta Delta$") + 
    scale_colour_manual(values = c("black", "grey")) +
    scale_linetype_manual(values = c(1, 1)) +
    theme_bw() +
    theme(axis.line = element_line(color='black'),
          plot.background = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.border = element_blank()) +
    theme(legend.position="best") )

# save plot
save_tikz_plot(p, width = 5, filename = "./tex/many-K.tex")

list.files('./data/results/joint-predictive', pattern = "jointpred*", full.names = TRUE) %>% 
  lapply(read_csv) %>% 
  bind_rows %>%
  write_csv('data/results/joint-predictive/jointpred_all.csv')

list.files('./data/results/many-irrelevant', pattern = "many_models*", full.names = TRUE) %>% 
  lapply(read_csv) %>% 
  bind_rows %>%
  write_csv('data/results/many-irrelevant/many_models_all.csv')

list.files('./data/results/all-irrelevant', pattern = "all_irrelevant*", full.names = TRUE) %>% 
  lapply(read_csv) %>% 
  bind_rows %>%
  write_csv('data/results/all-irrelevant/all_irrelevant_all.csv')

remotes::install_github("stan-dev/cmdstanr", Ncpus = 16)
install.packages(c("loo", "posterior", "tidyverse", "scoringRules"), Ncpus = 16)
remotes::install_github("TeemuSailynoja/bayesflow", Ncpus = 16)


## All irrelevant
## --------------

# read in results
results <- read.csv(file = "data/results/all_irrelevant_all.csv")

# plot results
plotting_df <- results %>% group_by(K, iter) %>%
  mutate(best_loo = max(loo_elpd_diff)) %>%
  group_by(K) %>%
  mutate(mean_best_loo = mean(best_loo)) %>%
  select(c("K", "mean_best_loo")) %>%
  unique()
ann_text_1 <- data.frame(K = 15, y = 3, lab = "Text")
ann_text_2 <- data.frame(K = 55, y = 8.5, lab = "Text")

# define order statistic heuristic data
order_stat_heuristic <- function(k, c = 2.8) {
  qnorm(p = 1 - 1 / (k * 2), mean = 0, sd = c)
}

# define the observed data
x <- plotting_df$K
y <- plotting_df$mean_best_loo
plot(y)

# define loss function
optim_fun <- function(par, y, x){
  # extract parameter
  c <- par[1]
  # compute loss
  y_pred <- x |> purrr::map_dbl(\(k) order_stat_heuristic(k, c = c))
  return(sum(y - y_pred)^2)
}

# fine tune the heuristic
optim_par <- optim(fn = optim_fun,
                 par = c(c = 1),
                 y = y,
                 x = x)
( c_optim <- optim_par$par )

# define plotting data
order_stat_df <- data.frame(K = x, 
                            y = x |> purrr::map_dbl(\(k) order_stat_heuristic(k, c = c_optim)))
(p <- ggplot() + 
    geom_line(data = order_stat_df, aes(x = K, y = y), 
              linetype = "dashed", colour = "grey") +
    geom_point(data = data.frame(x = x, y = y), aes(x = x, y = y)) +
    geom_label(data = ann_text_1, aes(x = K, y = y), label = "Eq. 21", 
               colour = "grey", size = 3) +
    geom_label(data = ann_text_2, aes(x = K, y = y), label = "Empirical observation", 
               colour = "black", size = 3) +
    #scale_x_continuous(trans='log2') +
    ylab("Best LOO-CV elpd difference") + 
    xlab("$K$") + 
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank()) )

# save plot
save_tikz_plot(p, width = 5, filename = "./tex/high-risk.tex")

