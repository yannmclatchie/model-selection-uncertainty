library(ggplot2)
library(geomtextpath)
library(dplyr)
library(bayesflow)

## Many irrelevant
## ---------------

# read in results
results <- read.csv(file = "data/old_results/many_models_all.csv")

# plot results
plotting_df <- results %>% group_by(K, beta, iter) %>%
  mutate(best_loo = max(loo_elpd_diff)) %>%
  group_by(K, beta) %>%
  mutate(mean_best_loo = median(best_loo)) %>%
  left_join(results %>% 
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
ann_text_1 <- data.frame(beta = 0.01,value = 3.5, lab = "Text",variable="loo",
                       K = factor(100,levels = c("2","10","100")))
ann_text_2 <- data.frame(beta = 0.01,value = 0.5,lab = "Text",variable="true",
                         K = factor(100,levels = c("2","10","100")))
scaleFUN <- function(x) sprintf("%.2f", x)

( p <- plotting_df %>% 
    ggplot(aes(x = beta, y = value, colour = variable)) + 
    geom_point() +
    geom_smooth(method = "gam", 
                formula = y ~ s(x, k = 6, bs = "cs", m = 1), 
                se = FALSE) +
    facet_wrap(vars(K), labeller = labeller(K = 
                                              c("2" = "$K = 2$",
                                                "10" = "$K = 10$",
                                                "100" = "$K = 100$")
    )) +
    geom_label(data = ann_text_1,label = "Best LOO-CV",size=3) +
    geom_label(data = ann_text_2,label = "True model",size=3,colour="grey") +
    scale_y_continuous(trans='pseudo_log') +
    scale_x_continuous(trans='log10', #labels=scaleFUN,
                       labels = function(x) ifelse(x == 0, "0", x)) +
    ylab("elpd difference") + 
    xlab("$beta Delta$") + 
    scale_colour_manual(values = c("black", "grey")) +
    scale_linetype_manual(values = c(1, 1)) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          panel.background = element_blank(),
          #panel.grid.major.y = element_line( size=.05, color="grey" )
          ) + 
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


devtools::load_all("../posterior")

# read in results
results <- read.csv(file = "data/all_irrelevant_all.csv") |>
  arrange(K, iter)
results

# define order statistic heuristic data
order_stat_heuristic <- function(k, c = 2.8) {
  qnorm(p = 1 - 1 / (k * 2), mean = 0, sd = c)
}

# plot results
plotting_df <- results |>
  group_by(K, iter) |>
  mutate(best_loo = max(loo_elpd_diff),
         mean_comb = mean(loo_elpd_diff),
         var_comb = sum(loo_elpd_diff_sd^2 + (loo_elpd_diff - mean_comb)^2) / K,
         sd_comb = sqrt(var_comb),
         diff_median = median(loo_elpd_diff),
         elpd_loo_diff_trunc = case_when(loo_elpd_diff >= diff_median ~ loo_elpd_diff - diff_median,
                                         loo_elpd_diff < diff_median ~ NA),
         n_models = sum(!is.na(elpd_loo_diff_trunc)),
         candidate_sd = sqrt(1 / n_models * sum(elpd_loo_diff_trunc^2, na.rm = TRUE)),
         order_stat = order_stat_heuristic(K, candidate_sd)) |>
  group_by(K) |>
  summarise(mean_best_loo = mean(best_loo),
            mean_order_stat = mean(order_stat)) |>
  ungroup() |>
  arrange(K)
plotting_df

pareto_order_stat <- function(n, sigma, k) {
  posterior:::qgeneralized_pareto(p = 1 - 1 / (n * 2), sigma = sigma, mu = 0) 
}

results |>
  dplyr::filter(K == 98) |>
  dplyr::filter(iter == 1) |>
  dplyr::summarise(khat = posterior::pareto_khat(loo_elpd_diff, r_eff = 1)$khat,
                   sigma = sd(loo_elpd_diff),
                   gauss_os = order_stat_heuristic(K, sigma),
                   pareto_os = pareto_order_stat(n = k, sigma = sigma, k = khat))

# smooth order statistics estimats
mono.spline <- scam::scam(mean_order_stat ~ s(K, k = 10, bs = "mpi", m = 1), 
                    data = plotting_df)
plotting_df$smoothed_order_stat <- mono.spline$fit

# define label data
ann_text_1 <- data.frame(K = 15, y = 2, lab = "Text")
ann_text_2 <- data.frame(K = 80, y = 1.5, lab = "Text")
p <- ggplot() + 
    geom_line(data = plotting_df, aes(x = K, y = smoothed_order_stat), 
              linetype = "dashed", colour = "grey") +
    geom_point(data = plotting_df, aes(x = K, y = mean_best_loo)) +
    geom_label(data = ann_text_1, aes(x = K, y = y), label = "Eq. 13", 
               colour = "grey", size = 3) +
    geom_label(data = ann_text_2, aes(x = K, y = y), label = "Empirical observation", 
               colour = "black", size = 3) +
    #scale_x_continuous(trans='log2') +
    #scale_y_continuous(trans='log2') +
    ylab("Best LOO-CV elpd difference") + 
    xlab("$K$") + 
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank())
p

# save plot
#save_tikz_plot(p, width = 5, filename = "./tex/high-risk.tex")

