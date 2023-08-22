library(ggplot2)
library(geomtextpath)
library(dplyr)
library(bayesflow)

## Many irrelevant
## ---------------

# read in results
results <- read.csv(file = "data/results/many_models_all.csv")

# pre-processing
results <- results |>
  mutate(test_elpd = test_elpd * n / n_test,
         baseline_test_elpd = baseline_test_elpd * n / n_test,
         loo_elpd_diff = loo_elpd - baseline_loo_elpd,
         test_elpd_diff = test_elpd - baseline_test_elpd) |>
  as_tibble()
results

# compute oracle model metrics
oracle_df <- results %>% 
  group_by(K, beta, model) %>%
  summarise(model_test_diffs = mean(test_elpd_diff)) %>%
  group_by(K, beta) %>%
  summarise(oracle_diff = max(max(model_test_diffs), 0))
oracle_df

# make plot
ann_text_1 <- data.frame(beta = 0.02, `value` = 6, lab = "Text",variable="selected-loo",
                         K = factor(100,levels = c("2","10","100")))
ann_text_2 <- data.frame(beta = 0.02, `value` = 1,lab = "Text",variable="true-test",
                         K = factor(100,levels = c("2","10","100")))
ann_text_3 <- data.frame(beta = 0.02, `value` = -1.6,lab = "Text",variable="selected-test",
                         K = factor(100,levels = c("2","10","100")))

p_oracle <- results %>% 
  slice_max(loo_elpd_diff, n = 1, by = c(K, beta, iter)) %>%
  group_by(K, beta) %>%
  summarise(mean_best_loo = mean(loo_elpd_diff),
            mean_best_test = mean(test_elpd_diff)) %>%
  left_join(results %>% 
            filter(model == 1) %>% # the true model
            group_by(K, beta) %>%
            summarise(mean_true_diff = mean(test_elpd_diff))) %>%
  #left_join(oracle_df) %>%
  #mutate(mean_best_loo = mean_best_loo - oracle_diff,
  #       mean_best_test = mean_best_test - oracle_diff,
  #       mean_true_diff = mean_true_diff - oracle_diff) |>
  #select(-oracle_diff) |>
  reshape2::melt(id = c("beta", "K")) %>%
  mutate_at(c("variable"), funs(recode(., `mean_best_loo`="selected-loo",
                                          `mean_best_test`="selected-test",
                                          `mean_true_diff`="true-test"))) %>%
  ggplot(aes(x = beta, y = value, colour = variable, linetype = variable)) +
  geom_hline(yintercept = 0, size = 0.2) +
  geom_point() +
  geom_smooth(method = "gam", 
              formula = y ~ s(x, k = 9, bs = "cs", m = 1), 
              se = F) +
  facet_wrap(vars(K), labeller = labeller(K = 
                                            c("2" = "$K = 2$",
                                              "10" = "$K = 10$",
                                              "100" = "$K = 100$")
  )) +
  geom_label(data = ann_text_1,label = "Selected\nLOO-CV",size=3) +
  geom_label(data = ann_text_2,label = "True test",size=3) +
  geom_label(data = ann_text_3,label = "Selected test",size=3) +
  #geom_label(data = ann_text_3,label = "Oracle test",size=3) +
  scale_y_continuous(trans='pseudo_log', breaks = c(-5, 0, 5, 100, 300)) +
  scale_x_continuous(trans='log10',
                     labels = function(x) ifelse(x == 0, "0", x)) +
  #ylab("elpd diff. to baseline model") + 
  ylab("elpd diff. to oracle model test elpd") + 
  xlab("$beta Delta$") + 
  scale_linetype_manual(values = c("solid", "solid","dashed", "solid")) +
  scale_colour_manual(values = c("red", "black", "grey", "blue")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.background = element_blank(),
        legend.position = "none"
  ) 
p_oracle
#save_tikz_plot(p_oracle, width = 5, filename = "./tex/many-K.tex")

## Data processing
## --------------

list.files('./data/results/many-irrelevant', pattern = "many_models*", full.names = TRUE) %>% 
  lapply(read_csv) %>% 
  bind_rows %>%
  write_csv('data/results/many-irrelevant/many_models_all.csv')

list.files('./data/results/all-irrelevant', pattern = "all_irrelevant*", full.names = TRUE) %>% 
  lapply(read_csv) %>% 
  bind_rows %>%
  write_csv('data/results/all-irrelevant/all_irrelevant_all.csv')

# check the Pareto k hats with posterior package
#remotes::install_github("stan-dev/cmdstanr", Ncpus = 16)
#install.packages(c("loo", "posterior", "tidyverse", "scoringRules"), Ncpus = 16)
#remotes::install_github("TeemuSailynoja/bayesflow", Ncpus = 16)
#devtools::load_all("../posterior")

## All irrelevant
## --------------

# read in results and pre-process
results <- read.csv(file = "data/results/all_irrelevant_all.csv") |>
  mutate(test_elpd = test_elpd * n / n_test,
         loo_elpd_diff = loo_elpd - baseline_loo_elpd)
results

# compute the mean best observed difference over iterations
mean_best_df <- results |>
  slice_max(loo_elpd_diff, n = 1, by = c(K, iter)) |>
  group_by(K) |>
  summarise(mean_best_diff = mean(loo_elpd_diff)) 

# define order statistic heuristic data
order_stat_heuristic <- function(k, c = 2.8) {
  qnorm(p = 1 - 1 / (k * 2), mean = 0, sd = c)
}

# compute the mean order statistic over iterations
mean_order_stat_df <- results |>
  group_by(K, iter) |>
  summarise(diff_median = median(loo_elpd_diff),
            elpd_loo_diff_trunc = case_when(loo_elpd_diff >= diff_median ~ loo_elpd_diff - diff_median,
                                            loo_elpd_diff < diff_median ~ NA),
            n_models = sum(!is.na(elpd_loo_diff_trunc)),
            candidate_sd = sqrt(1 / n_models * sum(elpd_loo_diff_trunc^2, na.rm = TRUE)),
            order_stat = order_stat_heuristic(K, candidate_sd)) |>
  group_by(K) |>
  summarise(mean_order_stat = mean(order_stat)) 

# bind results
plotting_df <- inner_join(mean_best_df, mean_order_stat_df)

# smooth order statistics estimates
mono.spline <- scam::scam(mean_order_stat ~ s(K, k = 10, bs = "mpi", m = 1), 
                          data = plotting_df)
plotting_df$smoothed_order_stat <- mono.spline$fit


# define label data
ann_text_1 <- data.frame(K = 15, y = 2.5, lab = "Text")
ann_text_2 <- data.frame(K = 60, y = 1.5, lab = "Text")
p <- ggplot() + 
    geom_line(data = plotting_df, aes(x = K, y = smoothed_order_stat), 
              linetype = "dashed", colour = "grey") +
    geom_point(data = plotting_df, aes(x = K, y = mean_best_diff)) +
    geom_label(data = ann_text_1, aes(x = K, y = y), label = "Equation 18", 
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
save_tikz_plot(p, width = 5, filename = "./tex/high-risk.tex")

