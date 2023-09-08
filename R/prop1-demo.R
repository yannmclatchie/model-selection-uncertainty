# load packages 
if(!requireNamespace("pacman"))install.packages("pacman")
pacman::p_load(purrr, brms, Matrix, dplyr, RColorBrewer, patchwork,
               loo, tidyr, multiverse, ggdist, ggplot2, stringr, 
               distributional, posterior, cmdstanr, haven, geomtextpath)
library(bayesflow)

# set # of cores 
nc <- parallel::detectCores() - 1

# load data 
path <- "data/datasets/eeg/pnas_povreductioneeg.dta"
data_eeg <- read_dta(path)

# data without NAs in absalpha
withoutNA <- !is.na(data_eeg$absalpha)
data_eeg <- data_eeg[withoutNA,]

# helper functions (from Andrew Gelman's code)
mean_impute <- function(a) ifelse(is.na(a), mean(a[!is.na(a)]), a)
standardize <- function(a) (a - mean(a))/(2*sd(a))

# create multiverse
# 1. pre_score definition
# 2. include different variables
# 3. normal model vs. log-normal model

# initialise multiverse
M_eeg = multiverse()

# define prior
prior = c(
  prior(horseshoe(df = 1, scale_slab = 1, df_slab = 100, par_ratio = 0.5), class = 'b')
)

inside(M_eeg, {
  
  # preprocessing (from Andrew Gelman's example on blog)
  df <- data_eeg %>%
    mutate(
      girl = as.numeric(cfemalea0),
      birthweight = as.numeric(cweightlba0_r),
      gestage = as.numeric(cgestagewksa0_r),
      momedu = mean_impute(as.numeric(momeduyrsa0)),
      income = mean_impute(as.numeric(hhrevisedincomea0)),
      white = as.numeric(mracea0) == 1,
      black = as.numeric(mracea0) == 2,
      momhealth = as.numeric(mgoodhealtha0),
      smoking = as.numeric(mcigduringavgwka0_r),
      drinking = as.numeric(malcduringavgwka0_r)
    )
  
  # 1. different pre_score definitions - pre_score1 from Andrews code
  df <- df %>%
    mutate(pre_score = branch(pre_score_calculation,"prescore_h" ~ standardize(birthweight) + standardize(gestage) +
                                standardize(momhealth) - standardize(smoking) - standardize(drinking)
    ))
  
  # 2. include different variables: treat, pre_score, girl 
  # 3. normal vs. log-normal model 
  mod_eeg <- brm(absalpha ~ 
                   branch(formula, 
                          "eq_1" ~ treat,
                          "eq_2" ~ treat + pre_score,
                          "eq_3" ~ treat + pre_score + income,
                          "eq_4" ~ treat + pre_score + birthweight,
                          "eq_5" ~ treat + pre_score + white + black
                   ),
                 data = df, 
                 prior = prior,
                 family = 
                   branch(family, 
                          "normal" ~ gaussian(), 
                          "lognormal" ~ lognormal()),
                 cores = nc,
                 backend = 'cmdstanr',
                 silent = 2,
                 refresh = 0)
})

inside(M_eeg, {
  # compute the model LOO-CV elpd
  mod_loo <- loo(mod_eeg)
  elpd <- mod_loo$estimates["elpd_loo", "Estimate"]
  se <- mod_loo$estimates["elpd_loo", "SE"]
  
  # extract the pointwise results
  pointwise_loo <- mod_loo$pointwise[,"elpd_loo"]
})

# execute the multiverse analysis
execute_multiverse(M_eeg)

loo_diff = function(loo1, loo2) {
  # compute the LOO-CV elpd difference between two LOO objecects
  if(is.null(loo1)) return(NA)
  if(is.null(loo2)) return(NA)
  stopifnot(nrow(loo1$pointwise) == nrow(loo2$pointwise))
  
  n = nrow(loo1$pointwise)
  diff_se = sd(loo1$pointwise[, 'elpd_loo'] - loo2$pointwise[, 'elpd_loo']) * sqrt(n)
  diff = sum(loo1$pointwise[, 'elpd_loo'] - loo2$pointwise[, 'elpd_loo'])
  return(list(elpd_loo_diff = diff, elpd_loo_diff_se = diff_se))
}

compute_elpd_diff = function(M_tibble) {
  # identify the best universe
  which_best <- M_tibble |>
    mutate( 
      elpd = map(.results, "elpd" ),
      se = map(.results, "se" )
    ) |>
    select( .universe, elpd, se ) |>
    unnest(elpd, se) |>
    filter(elpd == min(elpd)) |> # use WORST model as baseline
    pull(.universe)
  print(which_best)
  loo_df <- M_tibble |>
    mutate( 
      mod_loo = map(.results, "mod_loo" )
    ) |>
    select( .universe, mod_loo )
  best_loo <- loo_df |>
    filter(.universe == which_best) |>
    pull(mod_loo)
  best_loo <- best_loo[[1]]
  
  # compute the elpd differences
  elpd_diffs <- loo_df$mod_loo |> 
    map(\(mod_loo) loo_diff(mod_loo, best_loo))
  diff_df <- do.call(rbind, elpd_diffs) |>
    as_tibble() |>
    unnest(cols = c(elpd_loo_diff, elpd_loo_diff_se)) |>
    mutate(.universe = loo_df$.universe)
  
  # bind dataframes
  results <- diff_df |> 
    filter(elpd_loo_diff_se > 0)
  return(list(df = results, which_best = which_best))
}

# compute the summary statistics from the 
multiverse::expand(M_eeg)
M_tibble <- multiverse::expand(M_eeg)
results <- compute_elpd_diff(M_tibble)
diff_df <- results$df |> mutate(.data_name = "EEG cash aid")
diff_df

# concatenate data into plotting data
plot_df <- diff_df |> #rbind(diff_df, hier_df) |>
  add_row(elpd_loo_diff = 0, elpd_loo_diff_se = 0, 
          .universe = results$which_best, 
          .data_name = "EEG cash aid") 

# define order statistic heuristic data
order_stat_heuristic <- function(k, c = 2.8) {
  qnorm(p = 1 - 1 / (k * 2), mean = 0, sd = c)
}

# order statistics
order_stat_line_df <- plot_df |>
  summarise(diff_median = median(elpd_loo_diff),
            elpd_loo_diff_trunc = case_when(elpd_loo_diff >= diff_median ~ elpd_loo_diff - diff_median,
                                            elpd_loo_diff < diff_median ~ NA),
            n_models = sum(!is.na(elpd_loo_diff_trunc)),
            K = sum(!is.na(elpd_loo_diff)),
            candidate_sd = sqrt(1 / n_models * sum(elpd_loo_diff_trunc^2, na.rm = TRUE)),
            order_stat = order_stat_heuristic(K, candidate_sd)) |>
  select(order_stat) |>
  unique() |>
  mutate(.data_name = "EEG cash aid", label = "S(K)o K")
order_stat_line_df

# make plot
data_label_names = c(`0`='$rho = 0$', `0.5`='$rho = 0.5$', `0.9`='$rho = 0.9$')
p <- plot_df |>
  mutate(.universe = match(.universe, unique(.universe))) |>
  ggplot(aes(.universe, elpd_loo_diff)) +
  geom_hline(yintercept = 0, colour = "black", linetype = "solid", size = 0.2) +
  geom_pointrange(aes(ymin = elpd_loo_diff - elpd_loo_diff_se, 
                      ymax = elpd_loo_diff + elpd_loo_diff_se), 
                  shape = 1) +
  geom_labelhline(data = order_stat_line_df, 
                  aes(yintercept = order_stat, label = label), 
                  colour = "red", 
                  linetype = "dashed",
                  hjust = 0.1) +
  xlab("Model") +
  ylab("Delta elpdPlainHat") +
  scale_x_continuous(breaks = 1:10) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.background = element_blank(),
        legend.position="none") 
p
#save_tikz_plot(p, width = 4, filename = "./tex/prop1-demo.tex")
