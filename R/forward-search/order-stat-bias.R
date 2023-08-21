library(dplyr)
library(ggplot2)
library(bayesflow)

# load utilities
setwd("~/Desktop/model-selection-uncertainty")
source('R/forward-search/utils.R')
source('R/forward-search/forward_search.R')

# load the data
fs_r2d2 <- read.csv("data/results/fs.csv") 
fs_normal <- read.csv("data/results/fs_normal.csv") 
fs <- fs_r2d2 |> 
  rbind(fs_normal)
dim(fs)
candidates_r2d2 <- read.csv("data/results/fs_candidates.csv") 
candidates_normal <- read.csv("data/results/fs_normal_candidates.csv") |>
  mutate(name = "Normal prior")
candidates <- candidates_r2d2 |>
  rbind(candidates_normal)
dim(candidates)
n_pars <- 100

test_iter <- candidates |>
  filter(name == "r2d2", iter == 1, n == 100, rho == 0, size == 20)
  
norm_approx <- test_iter |>
  summarise(mu = mean(elpd_loo_diff),
            sd = sd(elpd_loo_diff)) 
norm_approx
half_norm_approx <- test_iter |>
  summarise(diff_median = median(elpd_loo_diff),
            elpd_loo_diff_trunc = case_when(elpd_loo_diff >= diff_median ~ elpd_loo_diff - diff_median,
                                            elpd_loo_diff < diff_median ~ NA),
            n_models = sum(!is.na(elpd_loo_diff_trunc)),
            half_sd = sqrt(1 / n_models * sum(elpd_loo_diff_trunc^2, na.rm = TRUE)),
            half_mu = diff_median) |>
  unique()
half_norm_approx

p <- test_iter |>
  ggplot(aes(elpd_loo_diff)) +
  geom_histogram(aes(y = ..density..), colour = "grey", fill = "grey") +
  geom_function(fun = dnorm, 
                args = list(mean = norm_approx$mu[1], 
                            sd = norm_approx$sd[1]), colour = "black") +
  geom_function(fun = dnorm, 
                args = list(mean = half_norm_approx$half_mu[1], 
                            sd = half_norm_approx$half_sd[1]), colour = "red") +
  xlab("Delta elpdHatPlain") +
  ylab("Density") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank())
p


# compute the combined standard deviation across candidates
fs_with_candidates <- candidates |> 
    group_by(name, iter, n, rho, size) |>
    summarise(diff_median = median(elpd_loo_diff),
              elpd_loo_diff_trunc = case_when(elpd_loo_diff >= diff_median ~ elpd_loo_diff - diff_median,
                                              elpd_loo_diff < diff_median ~ NA),
              n_models = sum(!is.na(elpd_loo_diff_trunc)),
              candidate_sd = sqrt(1 / n_models * sum(elpd_loo_diff_trunc^2, na.rm = TRUE))) |>
    ungroup() |> 
    select(name, iter, n, rho, size, diff_median, candidate_sd) |>
    unique() |>
    right_join(fs, by = c("name", "iter", "n", "rho", "size"))
fs_with_candidates

# compute the base loo-cv elpd
( fs_with_base <- fs_with_candidates |> 
    group_by(name, iter, n, rho) |>
    mutate(base_elpd_loo = first(elpd_loo)) |>
    ungroup() )

# define order statistic heuristic data
order_stat_heuristic <- function(k, mu = 0, sd = 1) {
  qnorm(p = 1 - 1 / (k * 2), mean = mu, sd = sd)
}

# compute the estimated bias
( fs_with_bias <- fs_with_base |> 
    #filter(size > 0) |>
    group_by(name, iter, n, rho) |>
    mutate(estimated_bias = order_stat_heuristic(k = n_pars + 1 - size, 
                                                 mu = 0,
                                                 sd = candidate_sd),
           corrected_diff = case_when(abs(elpd_loo_diff) > estimated_bias ~ elpd_loo_diff,
                                      abs(elpd_loo_diff) <= estimated_bias ~ elpd_loo_diff - 1.5 * estimated_bias)) |>
    ungroup() )

# correct the bias with the estimate
( fs_corrected <- fs_with_bias |>
    group_by(name, n, iter, rho) |>
    mutate(corrected_diff = ifelse(is.na(corrected_diff), 0, corrected_diff),
           corrected_loo = base_elpd_loo + cumsum(corrected_diff)) |>
    ungroup() )

# compute mlpds
( fs_mlpds <- fs_corrected |>
    group_by(name, iter, n, rho) |>
    mutate(mlpd_loo = (elpd_loo - elpd_loo_ref) / n, 
           mlpd_test = (elpd_test - elpd_test_ref) / n_test,
           mlpd_corrected = (corrected_loo - elpd_loo_ref) / n) |>
    ungroup() )

# compute the means over iterations
( means <-  fs_mlpds |>
  group_by(name, size, n, rho) |>
  summarize(mean_mlpd_loo = mean(mlpd_loo),
            mean_mlpd_corrected = mean(mlpd_corrected),
            mean_mlpd_test = mean(mlpd_test)) |>
  ungroup() )

# identify the point of maximal over-fitting
( over_fitting_df <- means |>
  group_by(name, n, rho) |>
  filter(mean_mlpd_loo == max(mean_mlpd_loo)) |>
  mutate(overfit_size = size) |>
  select(name, n, rho, overfit_size) |>
  ungroup() )

# remove heuristic calculation beyond point of maximal over-fitting
fs_mlpds <- fs_mlpds |> left_join(over_fitting_df) |>
  mutate(mlpd_corrected = ifelse(size > overfit_size, NA, mlpd_corrected))

# compute the means over iterations with new NA values
( means <-  fs_mlpds |>
    group_by(name, size, n, rho) |>
    summarize(mean_mlpd_loo = mean(mlpd_loo),
              mean_mlpd_corrected = mean(mlpd_corrected),
              mean_mlpd_test = mean(mlpd_test)) |>
    ungroup() )

# identify point of heuristic selection
order_stat_df <-  means |>
  group_by(name, n, rho) |>
  filter(mean_mlpd_corrected == max(mean_mlpd_corrected, na.rm = T)) |>
  mutate(overfit_size = size) |>
  select(name, n, rho, overfit_size) |>
  ungroup() 
order_stat_df

# plot the results
prior_to_plot <- "R2D2 prior" # "Normal prior"
prior_means <- means |> filter(name == prior_to_plot)
over_fitting_df_means <- over_fitting_df |> filter(name == prior_to_plot)
order_stat_df_means <- order_stat_df |> filter(name == prior_to_plot)
n_label_names = c(`100`='$n = 100$', `200`='$n = 200$', `400`='$n = 400$')
rho_label_names = c(`0`='$rho = 0$', `0.5`='$rho = 0.5$', `0.9`='$rho = 0.9$')
p <- fs_mlpds |>
  filter(name == prior_to_plot,
         rho != 0.5) |>
  ggplot() +
  geom_hline(yintercept = 0, size = 0.2) + 
  geom_line(aes(x = size, y = mlpd_loo, group = iter, color = 'Train (MLPD LOO)'), 
            size = 0.1, 
            alpha = 0.15) +
  geom_line(aes(x = size, y = mlpd_test, group = iter, color = 'Test (MLPD)'), 
            size = 0.1, 
            alpha = 0.15) +
  geom_line(aes(x = size, y = mlpd_corrected, group = iter, color = 'Corrected (MLPD LOO)'), 
            size = 0.1, 
            alpha = 0.15) +
  geom_line(data = prior_means |> filter(rho != 0.5),
            aes(x = size, y = mean_mlpd_loo, color = 'Train (MLPD LOO)'),
            size=0.9) +
  geom_line(data = prior_means |> filter(rho != 0.5),
            aes(x = size, y = mean_mlpd_test, color = 'Test (MLPD)'),
            size=0.9) +
  geom_line(data = prior_means |> filter(rho != 0.5),
            aes(x = size, y = mean_mlpd_corrected, color = 'Corrected (MLPD LOO)'),
            size=0.9) +
  geom_vline(data = order_stat_df_means |> filter(rho != 0.5), 
             aes(xintercept = overfit_size), 
             size = 0.5, colour = 'black', linetype = 'dashed') +
  geom_vline(data = over_fitting_df_means|> filter(rho != 0.5),
             aes(xintercept = overfit_size), 
             size = 0.5, colour = 'red', linetype = 'dashed') +
  scale_colour_manual("", 
                      breaks = c("Train (MLPD LOO)", "Test (MLPD)", "Corrected (MLPD LOO)"),
                      values = c("red", "grey", "black")) +
  facet_grid(
    rows = vars(n), 
    cols = vars(rho), 
    labeller = labeller(n = as_labeller(n_label_names), 
                        rho = as_labeller(rho_label_names))) +
  ylim(-0.5, 0.5) +
  xlab("Model size") +
  ylab("mlpd") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.background = element_blank(),
        legend.position="none") 
p
#save_tikz_plot(p, width = 5, filename = "./tex/gaussian-forward.tex")
#save_tikz_plot(p, width = 5, filename = "./tex/r2d2-forward.tex")


##
## Using different asymptotic correlation assumptions
##

# compute the estimated bias
( fs_with_bias_multi <- fs_with_base |> 
    #filter(size > 0) |>
    group_by(name, iter, n, rho) |>
    mutate(estimated_bias = order_stat_heuristic(k = n_pars + 1 - size, 
                                                 mu = 0,
                                                 sd = candidate_sd),
           corrected_diff_1 = case_when(abs(elpd_loo_diff) > estimated_bias ~ elpd_loo_diff,
                                      abs(elpd_loo_diff) <= estimated_bias ~ elpd_loo_diff - 1 * estimated_bias),
           corrected_diff_15 = case_when(abs(elpd_loo_diff) > estimated_bias ~ elpd_loo_diff,
                                        abs(elpd_loo_diff) <= estimated_bias ~ elpd_loo_diff - 1.5 * estimated_bias),
           corrected_diff_2 = case_when(abs(elpd_loo_diff) > estimated_bias ~ elpd_loo_diff,
                                        abs(elpd_loo_diff) <= estimated_bias ~ elpd_loo_diff - 2 * estimated_bias)) |>
    ungroup() )

# correct the bias with the estimate
( fs_corrected_multi <- fs_with_bias_multi |>
    group_by(name, n, iter, rho) |>
    mutate(corrected_diff_1 = ifelse(is.na(corrected_diff_1), 0, corrected_diff_1),
           corrected_diff_15 = ifelse(is.na(corrected_diff_15), 0, corrected_diff_15),
           corrected_diff_2 = ifelse(is.na(corrected_diff_2), 0, corrected_diff_2),
           corrected_loo_1 = base_elpd_loo + cumsum(corrected_diff_1),
           corrected_loo_15 = base_elpd_loo + cumsum(corrected_diff_15),
           corrected_loo_2 = base_elpd_loo + cumsum(corrected_diff_2)) |>
    ungroup() )

# compute mlpds
( fs_mlpds_multi <- fs_corrected_multi |>
    group_by(name, iter, n, rho) |>
    mutate(mlpd_loo = (elpd_loo - elpd_loo_ref) / n, 
           mlpd_test = (elpd_test - elpd_test_ref) / n_test,
           mlpd_corrected_1 = (corrected_loo_1 - elpd_loo_ref) / n,
           mlpd_corrected_15 = (corrected_loo_15 - elpd_loo_ref) / n,
           mlpd_corrected_2 = (corrected_loo_2 - elpd_loo_ref) / n) |>
    ungroup() )

# compute the means over iterations
( means_multi <-  fs_mlpds_multi |>
    group_by(name, size, n, rho) |>
    summarize(mean_mlpd_loo = mean(mlpd_loo),
              mean_mlpd_corrected_1 = mean(mlpd_corrected_1),
              mean_mlpd_corrected_15 = mean(mlpd_corrected_15),
              mean_mlpd_corrected_2 = mean(mlpd_corrected_2),
              mean_mlpd_test = mean(mlpd_test)) |>
    ungroup() )

# identify the point of maximal over-fitting
( over_fitting_df_multi <- means_multi |>
    group_by(name, n, rho) |>
    filter(mean_mlpd_loo == max(mean_mlpd_loo)) |>
    mutate(overfit_size = size) |>
    select(name, n, rho, overfit_size) |>
    ungroup() )

# remove heuristic calculation beyond point of maximal over-fitting
fs_mlpds_multi <- fs_mlpds_multi |> left_join(over_fitting_df) |>
  mutate(mlpd_corrected_1 = ifelse(size > overfit_size, NA, mlpd_corrected_1),
         mlpd_corrected_15 = ifelse(size > overfit_size, NA, mlpd_corrected_15),
         mlpd_corrected_2 = ifelse(size > overfit_size, NA, mlpd_corrected_2))

# compute the means over iterations with new NA values
( means_multi <-  fs_mlpds_multi |>
    group_by(name, size, n, rho) |>
    summarize(mean_mlpd_loo = mean(mlpd_loo),
              mean_mlpd_corrected_1 = mean(mlpd_corrected_1),
              mean_mlpd_corrected_15 = mean(mlpd_corrected_15),
              mean_mlpd_corrected_2 = mean(mlpd_corrected_2),
              mean_mlpd_test = mean(mlpd_test)) |>
    ungroup() )

# plot the results
prior_to_plot <- "Normal prior" # "Normal prior"
prior_means <- means |> filter(name == prior_to_plot)
over_fitting_df_means <- over_fitting_df |> filter(name == prior_to_plot)
order_stat_df_means <- order_stat_df |> filter(name == prior_to_plot)
n_label_names = c(`100`='$n = 100$', `200`='$n = 200$', `400`='$n = 400$')
rho_label_names = c(`0`='$rho = 0$', `0.5`='$rho = 0.5$', `0.9`='$rho = 0.9$')
p_alt <- means_multi |>
  filter(name == prior_to_plot,
         rho != 0.5) |>
  ggplot() +
  geom_hline(yintercept = 0, size = 0.2) + 
  geom_line(aes(x = size, y = mean_mlpd_loo), colour = "red",
            size=0.9, alpha = 0.2) +
  geom_line(aes(x = size, y = mean_mlpd_test), color = 'blue',
            size=0.9, alpha = 0.2) +
  geom_line(aes(x = size, y = mean_mlpd_corrected_1),
            size=1.2, alpha = 0.4, colour = "black") +
  geom_line(aes(x = size, y = mean_mlpd_corrected_15),
            size=1.2, alpha = 0.6, colour = "black") +
  geom_line(aes(x = size, y = mean_mlpd_corrected_2),
            size=1.2, alpha = 1, colour = "black") +
  #geom_vline(data = order_stat_df_means |> filter(rho != 0.5), 
  #           aes(xintercept = overfit_size), 
  #           size = 0.5, colour = 'black', linetype = 'dashed') +
  #geom_vline(data = over_fitting_df_means|> filter(rho != 0.5),
  #           aes(xintercept = overfit_size), 
  #           size = 0.5, colour = 'red', linetype = 'dashed') +
  facet_grid(
    rows = vars(n), 
    cols = vars(rho), 
    labeller = labeller(n = as_labeller(n_label_names), 
                        rho = as_labeller(rho_label_names))) +
  ylim(-0.5, 0.5) +
  xlim(0, 50) +
  xlab("Model size") +
  ylab("mlpd") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.background = element_blank(),
        legend.position="none") 
p_alt
#save_tikz_plot(p_alt, width = 5, filename = "./tex/alt-gaussian-forward.tex")
#save_tikz_plot(p_alt, width = 5, filename = "./tex/alt-r2d2-forward.tex")
