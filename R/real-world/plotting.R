library(dplyr)
library(purrr)
library(ggplot2)
library(bayesflow)
library(tidyr)

## Data reading

# define data reading utility file
read_data = function(what, path) {
  stopifnot(what %in% c('projpred', 'forward_search', 
                        'forward_search_normal', 'candidates',
                        'candidates_normal'))
  files = dir(path)
  res = list()
  for (file in files) {
    d = readRDS(file.path(path, file))
    
    if(what == 'candidates' || what == 'candidates_normal') {
      candidates = d[[what]]
      fs = d[['forward_search']]
      n = unique(fs$n)[1]
      fold = unique(fs$fold)[1]
      data_name = unique(fs$data_name)[1]
      prior_name = unique(fs$prior_name)[1]
      candidates = map2(candidates, seq_along(candidates)-1, function(x, y) {
        x$size = y
        x$n = n
        x$fold = fold
        x$data_name = data_name
        x$prior_name = prior_name
        x
      }
      )
      candidates = do.call(rbind,candidates)
      res[[length(res)+1]] = candidates
    } else {
      res[[length(res)+1]] = d[[what]] 
    }
  }
  do.call(rbind, res) %>% ungroup()
}

# read in files
files <- list.files("data/results/real", full.names = TRUE)

# read R2D2 forward search results
fs_r2d2 <- files |> 
  map(readRDS) |> 
  map(rbind) |> 
  map(\(x) x[,"forward_search"]$forward_search) |> 
  bind_rows()
fs_r2d2

# read normal forward search results
fs_normal <- files |> 
  map(readRDS) |> 
  map(rbind) |> 
  map(\(x) x[,"forward_search_normal"]$forward_search_normal) |> 
  bind_rows() 
fs_normal

# concatenate R2D2 and normal forward search results
fs <- fs_r2d2 |> 
  rbind(fs_normal)
fs |> select(data_name) |> unique()

# read r2d2 candidates
r2d2_candidates <- read_data('candidates', 'data/results/real')

# read normal candidates
normal_candidates <- read_data('candidates_normal', 'data/results/real') |>
  mutate(prior_name = "Normal prior")

# concatenate the candidates
candidates <- r2d2_candidates |>
  rbind(normal_candidates)
candidates |> select(data_name) |> unique()

## Plotting

# compute the combined standard deviation across candidates
fs_with_candidates <- candidates |> 
  group_by(data_name, prior_name, fold, size) |>
  summarise(diff_median = median(elpd_loo_diff),
            elpd_loo_diff_trunc = case_when(elpd_loo_diff >= diff_median ~ elpd_loo_diff - diff_median,
                                            elpd_loo_diff < diff_median ~ NA),
            n_models = sum(!is.na(elpd_loo_diff_trunc)),
            candidate_sd = sqrt(1 / n_models * sum(elpd_loo_diff_trunc^2, na.rm = TRUE))) |>
  ungroup() |> 
  select(data_name, prior_name, fold, size, diff_median, candidate_sd) |>
  unique() |>
  right_join(fs, by = c("data_name", "prior_name", "fold", "size"))
fs_with_candidates

fs_with_candidates |> filter(prior_name == "Normal prior") |>
  select(fold, size, diff_median, candidate_sd, elpd_loo_diff)

# compute the base loo-cv elpd
( fs_with_base <- fs_with_candidates |> 
    group_by(data_name, prior_name, fold) |>
    mutate(base_elpd_loo = first(elpd_loo)) |>
    ungroup() )

# define order statistic heuristic data
order_stat_heuristic <- function(k, mu = 0, sd = 1) {
  qnorm(p = 1 - 1 / (k * 2), mean = mu, sd = sd)
}

# compute the estimated bias
fs_with_bias <- fs_with_base |> 
  mutate(n_pars = case_when(data_name == "sonar" ~ 60,
                            data_name == "ionosphere" ~ 32,
                            data_name == "colon" ~ 2000,
                            data_name == "crime" ~ 102,
                            data_name == "heart" ~ 44)) |>
  group_by(data_name, prior_name, fold, size) |>
    mutate(estimated_bias = order_stat_heuristic(k = n_pars + 1 - size, 
                                                 mu = 0,
                                                 sd = candidate_sd),
           corrected_diff = case_when(abs(elpd_loo_diff) > estimated_bias ~ elpd_loo_diff,
                                      abs(elpd_loo_diff) <= estimated_bias ~ diff_median)) |>
    ungroup()
fs_with_bias

fs_with_bias |> 
  filter(prior_name == "Normal prior") |>
  select(corrected_diff, estimated_bias, candidate_sd, elpd_loo_diff) |>
  print(n = 30)

# correct the bias with the estimate
fs_corrected <- fs_with_bias |>
  group_by(data_name, prior_name, fold) |>
    mutate(corrected_diff = ifelse(is.na(corrected_diff), 0, corrected_diff),
           corrected_loo = base_elpd_loo + cumsum(corrected_diff)) |>
    ungroup()
fs_corrected

fs_corrected |> 
  filter(prior_name == "Normal prior") |>
  select(base_elpd_loo, corrected_diff, corrected_loo, elpd_loo)

# compute mlpds
fs_mlpds <- fs_corrected |>
  group_by(data_name, prior_name, fold, size) |>
    mutate(mlpd_loo = (elpd_loo - elpd_loo_ref) / n, 
           mlpd_test = (elpd_test - elpd_test_ref) / n_test,
           mlpd_corrected = (corrected_loo - elpd_loo_ref) / n) |>
    ungroup()
fs_mlpds

# compute the means over folds
means <-  fs_mlpds |>
  group_by(data_name, prior_name, size) |>
  summarize(mean_mlpd_loo = mean(mlpd_loo),
            mean_mlpd_corrected = mean(mlpd_corrected),
            mean_mlpd_test = mean(mlpd_test)) |>
  ungroup()
means

# identify the point of maximal over-fitting
over_fitting_df <- means |>
  filter(data_name != 'colon') |>
  group_by(data_name, prior_name) |>
    filter(mean_mlpd_loo == max(mean_mlpd_loo)) |>
    mutate(overfit_size = size) |>
    select(data_name, prior_name, overfit_size) |>
    ungroup()
over_fitting_df

# remove heuristic calculation beyond point of maximal over-fitting
fs_mlpds <- fs_mlpds |> left_join(over_fitting_df) |>
  mutate(mlpd_corrected = ifelse(size > overfit_size, NA, mlpd_corrected))
fs_mlpds

# compute the means over iterations with new NA values
means <-  fs_mlpds |>
  filter(data_name != 'colon') |>
  group_by(data_name, prior_name, size) |>
    summarize(mean_mlpd_loo = mean(mlpd_loo),
              mean_mlpd_corrected = mean(mlpd_corrected),
              mean_mlpd_test = mean(mlpd_test)) |>
    ungroup()
means

# identify point of heuristic selection
order_stat_df <-  means |>
  filter(data_name != 'colon') |>
  group_by(data_name, prior_name) |>
  filter(mean_mlpd_corrected == max(mean_mlpd_corrected, na.rm = T)) |>
  mutate(overfit_size = size) |>
  select(data_name, prior_name, overfit_size) |>
  ungroup() 
order_stat_df

# make plot
label_names = c(`ionosphere`='Ionosphere', `heart`='Heart', `sonar`='Sonar', `crime`='Crime', `colon`='Colon')
prior_names = c(`R2D2 prior`='Sparsity-inducing prior', `Normal prior`='Normal prior')
p <- fs_mlpds |>
  filter(data_name != 'colon') |>
  ggplot() +
  geom_hline(yintercept = 0, size = 0.2) + 
  geom_line(aes(x = size, y = mlpd_loo, group = fold, color = 'Train (MLPD LOO)'), 
            size = 0.1, 
            alpha = 0.15) +
  geom_line(aes(x = size, y = mlpd_test, group = fold, color = 'Test (MLPD)'), 
            size = 0.1, 
            alpha = 0.15) +
  geom_line(aes(x = size, y = mlpd_corrected, group = fold, color = 'Corrected (MLPD LOO)'), 
            size = 0.1, 
            alpha = 0.15) +
  geom_line(data = means,
            aes(x = size, y = mean_mlpd_loo, color = 'Train (MLPD LOO)'),
            size=0.9) +
  geom_line(data = means,
            aes(x = size, y = mean_mlpd_test, color = 'Test (MLPD)'),
            size=0.9) +
  geom_line(data = means,
            aes(x = size, y = mean_mlpd_corrected, color = 'Corrected (MLPD LOO)'),
            size=0.9) +
  geom_vline(data = order_stat_df, aes(xintercept = overfit_size), 
             size = 0.5, colour = 'black', linetype = 'dashed') +
  geom_vline(data = over_fitting_df, aes(xintercept = overfit_size), 
             size = 0.5, colour = 'red', linetype = 'dashed') +
  scale_colour_manual("", 
                      breaks = c("Train (MLPD LOO)", "Test (MLPD)", "Corrected (MLPD LOO)"),
                      values = c("red", "grey", "black")) +
  facet_grid(
    rows = vars(data_name), 
    cols = vars(prior_name),
    scales = "free",
    labeller = labeller(data_name = as_labeller(label_names),
                        prior_name = as_labeller(prior_names))) +
  ylim(-0.4, 0.2) +
  xlim(0, 50) +
  xlab("Model size") +
  ylab("mlpd") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.background = element_blank(),
        legend.position="none") 
p
#bayesflow::save_tikz_plot(p, width = 5, height = 5 * 1.2, 
#                          filename = "./tex/real-world-forward.tex")

## Size selection

# read in projpred files
out_projpred <- files |> 
  map(readRDS) |> 
  map(rbind) |> 
  map(\(x) x[,"projpred"]$projpred) |> 
  bind_rows()
out_projpred

# compute the summary statistics of projpred sizes
projpred_sizes <- out_projpred |>
  mutate(elpd_loo = elpd_loo_se,
         elpd_loo_se = elpd_loo_diff_ref,
         elpd_loo_diff_ref = elpd_loo_diff_ref_se,
         elpd_loo_diff_ref_se = NA.x,
         elpd_test = elpd_test_se,
         elpd_test_se = elpd_test_diff_ref,
         elpd_test_diff_ref = elpd_test_diff_ref_se,
         elpd_test_diff_ref_se = NA.y,
         upper_ci = elpd_loo_diff_ref + elpd_loo_diff_ref_se) |>
  group_by(data_name, prior_name, fold) |>
  #summarise(size = first(which(upper_ci > 0))) |>
  summarise(size = first(which(elpd_loo_diff_ref > -4))) |>
  ungroup() 
projpred_sizes |> print(n = 30)

projpred_elpds <- out_projpred |> 
  mutate(elpd_test = elpd_test_se,
         mlpd_diff = (elpd_test - elpd_test_ref) / n_test)
projpred_elpds |> select(size, mlpd_diff) |> print(n = 30)

projpred_sizes_with_elpd <- projpred_sizes |>
  left_join(projpred_elpds, by = c("data_name", "prior_name", "fold", "size")) |>
  select(data_name, prior_name, size, mlpd_diff) |>
  pivot_longer(cols = c(size, mlpd_diff),
               names_to = "metric", values_to = "value")
projpred_sizes_with_elpd

projpred_metrics <- projpred_sizes_with_elpd |>
  group_by(data_name, prior_name, metric) |>
  summarise(mean_metric = mean(value, na.rm = T),
            lower_metric = quantile(value, probs = 0.05, na.rm = T),
            upper_metric = quantile(value, probs = 0.95, na.rm = T)) |>
  ungroup() |>
  mutate(heuristic = "Projpred", prior_name = "R2D2 prior")
projpred_metrics

# compute the bulge size
bulge_size <- fs |>
  group_by(data_name, prior_name, fold) |>
  filter(elpd_loo == max(elpd_loo)) |>
  select(data_name, prior_name, fold, size) |>
  ungroup() |>
  mutate(heuristic = "Bulge") |>
  left_join(fs, by = c("data_name", "prior_name", "fold", "size")) |>
  mutate(mlpd_diff = (elpd_test - elpd_test_ref) / n_test) |>
  select(data_name, prior_name, fold, size, mlpd_diff) |>
  pivot_longer(cols = c(size, mlpd_diff),
               names_to = "metric", values_to = "value")
bulge_size

bulge_metrics <- bulge_size |>
  group_by(data_name, prior_name, metric) |>
  summarise(mean_metric = mean(value, na.rm = T),
            lower_metric = quantile(value, probs = 0.05, na.rm = T),
            upper_metric = quantile(value, probs = 0.95, na.rm = T)) |>
  ungroup() |>
  mutate(heuristic = "Bulge")
bulge_metrics

# compute the size of smallest model indistinguishable from bulge
bulge_elpd <- fs |>
  group_by(data_name, prior_name, fold) |>
  filter(elpd_loo == max(elpd_loo)) |>
  select(data_name, prior_name, fold, size) |>
  ungroup() |>
  mutate(heuristic = "Bulge") |>
  left_join(fs, by = c("data_name", "prior_name", "fold", "size")) |>
  select(data_name, prior_name, fold, size, elpd_loo) |> 
  mutate(bulge_elpd = elpd_loo, bulge_size = size)
bulge_elpd

smallest_bulge_size <- fs |>
  mutate(upper_tail = elpd_loo + 2 * elpd_loo_se) |>
  left_join(bulge_elpd, by = c("data_name", "prior_name", "fold")) |>
  group_by(data_name, prior_name, fold) |>
  summarise(size = first(which(upper_tail >= bulge_elpd)) - 1) |>
  ungroup() |>
  mutate(heuristic = "$2s$") |>
  left_join(fs, by = c("data_name", "prior_name", "fold", "size")) |>
  mutate(mlpd_diff = (elpd_test - elpd_test_ref) / n_test) |>
  select(data_name, prior_name, fold, size, mlpd_diff) |>
  pivot_longer(cols = c(size, mlpd_diff),
               names_to = "metric", values_to = "value")
smallest_bulge_size

smallest_bulge_metrics <- smallest_bulge_size |>
  group_by(data_name, prior_name, metric) |>
  summarise(mean_metric = mean(value, na.rm = T),
            lower_metric = quantile(value, probs = 0.05, na.rm = T),
            upper_metric = quantile(value, probs = 0.95, na.rm = T)) |>
  ungroup() |>
  mutate(heuristic = "$2s$")
smallest_bulge_metrics

# compute heuristic size selection
order_stat_sizes <- fs_corrected |>
  group_by(data_name, prior_name, fold) |>
  filter(corrected_loo == max(corrected_loo, na.rm = T)) |>
  select(data_name, prior_name, fold, size) |>
  ungroup() |>
  mutate(heuristic = "Ours") |>
  left_join(fs, by = c("data_name", "prior_name", "fold", "size")) |>
  mutate(mlpd_diff = (elpd_test - elpd_test_ref) / n_test) |>
  select(data_name, prior_name, fold, size, mlpd_diff) |> 
  pivot_longer(cols = c(size, mlpd_diff),
               names_to = "metric", values_to = "value")
order_stat_sizes

order_stat_metrics <- order_stat_sizes |>
  group_by(data_name, prior_name, metric) |>
  summarise(mean_metric = mean(value, na.rm = T),
            lower_metric = quantile(value, probs = 0.05, na.rm = T),
            upper_metric = quantile(value, probs = 0.95, na.rm = T)) |>
  ungroup() |>
  mutate(heuristic = "Ours")
order_stat_metrics

# combine and plot
label_names = c(`ionosphere`='Ionosphere', `heart`='Heart', `sonar`='Sonar', `crime`='Crime', `colon`='Colon')
metric_names = c(`mlpd_diff`='$Delta mlpd$', `size`='Selected size')
p_sizes <- rbind(projpred_metrics, order_stat_metrics, smallest_bulge_metrics, bulge_metrics) |>
  filter(data_name != "colon") |>
  mutate(prior_name = case_when(prior_name == "Normal prior" ~ "Normal", 
                                prior_name == "R2D2 prior" ~ "R2D2")) |>
  ggplot(aes(y = mean_metric, x = heuristic, colour = prior_name)) +
  geom_hline(data = data.frame(xint = 0, 
                               prior_name = c('crime', 'heart', 'ionosphere', 'sonar'),
                               metric = 'mlpd_diff'), 
             aes(yintercept = xint), 
             size = 0.2) +
  geom_pointrange(aes(ymin = lower_metric, 
                      ymax = upper_metric),
                  size = 0.4,
                  shape = 1,
                  position = position_dodge(width = 0.4)) +
  facet_grid(vars(metric), vars(data_name), scales = "free",
             labeller = labeller(data_name = as_labeller(label_names),
                                 metric = as_labeller(metric_names))) +
  xlab("Heuristic") +
  ylab("") +
  scale_colour_manual(values = c("black", "red")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust=1),
        axis.title.y=element_blank())
p_sizes
bayesflow::save_tikz_plot(p_sizes, width = 5, filename = "./tex/real-world-sizes.tex")

## Appendix plots

# compute the 2 sigma metrics
fs_with_ci <- candidates |> 
  group_by(data_name, prior_name, fold, size) |>
  mutate(lower_2sigma_ci = elpd_loo_diff - 2 * elpd_loo_diff_se,
         lower_3sigma_ci = elpd_loo_diff - 3 * elpd_loo_diff_se) |>
  summarise(best_2sigma_ci = max(lower_2sigma_ci),
            best_3sigma_ci = max(lower_3sigma_ci)) |>
  ungroup() |>
  right_join(fs, by = c("data_name", "prior_name", "fold", "size"))
twosigma_sizes <- fs_with_ci |>
  group_by(data_name, prior_name, fold) |>
  summarise(size = first(which(best_2sigma_ci <= 0)) - 1) |>
  ungroup() |>
  mutate(heuristic = "$2s$") |>
  left_join(fs, by = c("data_name", "prior_name", "fold", "size")) |>
  mutate(mlpd_diff = (elpd_test - elpd_test_ref) / n_test) |>
  select(data_name, prior_name, fold, size, mlpd_diff) |>
  pivot_longer(cols = c(size, mlpd_diff),
               names_to = "metric", values_to = "value")
twosigma_metrics <- twosigma_sizes |>
  group_by(data_name, prior_name, metric) |>
  summarise(mean_metric = mean(value, na.rm = T),
            lower_metric = quantile(value, probs = 0.05, na.rm = T),
            upper_metric = quantile(value, probs = 0.95, na.rm = T)) |>
  ungroup() |>
  mutate(heuristic = "$2s$")
twosigma_metrics

# compute the 3 sigma metrics
threesigma_sizes <- fs_with_ci |>
  group_by(data_name, prior_name, fold) |>
  summarise(size = first(which(best_3sigma_ci <= 0)) - 1) |>
  ungroup() |>
  mutate(heuristic = "$3s$") |>
  left_join(fs, by = c("data_name", "prior_name", "fold", "size")) |>
  mutate(mlpd_diff = (elpd_test - elpd_test_ref) / n_test) |>
  select(data_name, prior_name, fold, size, mlpd_diff) |>
  pivot_longer(cols = c(size, mlpd_diff),
               names_to = "metric", values_to = "value")
threesigma_metrics <- threesigma_sizes |>
  group_by(data_name, prior_name, metric) |>
  summarise(mean_metric = mean(value, na.rm = T),
            lower_metric = quantile(value, probs = 0.05, na.rm = T),
            upper_metric = quantile(value, probs = 0.95, na.rm = T)) |>
  ungroup() |>
  mutate(heuristic = "$3s$")
threesigma_metrics

# make plot
label_names = c(`ionosphere`='Ionosphere', `heart`='Heart', `sonar`='Sonar', `crime`='Crime', `colon`='Colon')
metric_names = c(`mlpd_diff`='$Delta mlpd$', `size`='Selected size')
p_appendix <- rbind(order_stat_metrics, twosigma_metrics, 
                    threesigma_metrics) |>
  filter(data_name != 'colon') |>
  mutate(prior_name = case_when(prior_name == "Normal prior" ~ "Normal", 
                                prior_name == "R2D2 prior" ~ "R2D2")) |>
  ggplot(aes(y = mean_metric, x = heuristic, colour = prior_name)) +
  geom_hline(data = data.frame(xint = 0, 
                               prior_name = c('crime', 'heart', 'ionosphere', 'sonar'),
                               metric = 'mlpd_diff'), 
             aes(yintercept = xint), 
             size = 0.2) +
  geom_pointrange(aes(ymin = lower_metric, 
                      ymax = upper_metric),
                  size = 0.4,
                  shape = 1,
                  position = position_dodge(width = 0.4)) +
  facet_grid(vars(metric), vars(data_name), scales = "free",
             labeller = labeller(data_name = as_labeller(label_names),
                                 metric = as_labeller(metric_names))) +
  xlab("Heuristic") +
  ylab("") +
  scale_colour_manual(values = c("black", "red")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust=1),
        axis.title.y=element_blank())
p_appendix
bayesflow::save_tikz_plot(p_appendix, width = 5, filename = "./tex/incremental-sizes.tex")

