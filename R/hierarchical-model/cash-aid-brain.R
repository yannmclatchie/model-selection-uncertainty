# load packages 
if(!requireNamespace("pacman"))install.packages("pacman")
pacman::p_load(purrr, brms, Matrix, dplyr, RColorBrewer, patchwork,
               loo, tidyr, multiverse, ggdist, ggplot2, stringr, 
               distributional, posterior, cmdstanr, haven)
library(bayesflow)

# set # of cores 
nc <- parallel::detectCores() - 1

# load data 
setwd("~/Desktop/model-selection-uncertainty")
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
    mutate(pre_score = branch(pre_score_calculation,
                              "prescore_ag" ~ standardize(birthweight) + standardize(gestage) + 
                                standardize(momedu) + standardize(income) + white - black + 
                                standardize(momhealth) - standardize(smoking) - standardize(drinking)
                              #"pre_score2" ~ standardize(momedu) + standardize(income) + white - black,
                              #"prescore_h" ~ standardize(birthweight) + standardize(gestage) +
                              #standardize(momhealth) - standardize(smoking) - standardize(drinking)
    ))
  
  # 2. include different variables: treat, pre_score, girl 
  # 3. normal vs. log-normal model 
  mod_eeg <- brm(absalpha ~ 
                   branch(formula, 
                          "eq_1" ~ treat,
                          "eq_2" ~ treat + pre_score,
                          "eq_3" ~ treat + pre_score + girl,
                          "eq_4" ~ treat + girl + birthweight + gestage + 
                            momedu + income + white + black + momhealth + 
                            smoking + drinking
                   ),
                 data = df, 
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

compute_elpd_diff = function(M) {
  # identify the best universe
  which_best <- multiverse::expand(M) |>
    mutate( 
      elpd = map(.results, "elpd" ),
      se = map(.results, "se" )
    ) |>
    select( .universe, elpd, se ) |>
    unnest(elpd, se) |>
    filter(elpd == max(elpd)) |>
    pull(.universe)
  print(which_best)
  loo_df <- multiverse::expand(M) |>
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
    mutate(.universe = loo_df$.universe,
           .method = "raw")
  
  # compute the correlation between pointwise elpds
  elpd_mat <- multiverse::expand(M) |>
    mutate( 
      pointwise_loo = map(.results, "pointwise_loo" )
    ) |>
    select( .universe, pointwise_loo ) |>
    unnest(pointwise_loo) |>
    ungroup() |>
    pivot_wider(names_from = .universe, values_from = pointwise_loo) |>
    unnest() |>
    as.matrix()
  Rho <- cor(elpd_mat)
  
  # extract the best elpd by point estimate and its standard error
  se_best <- multiverse::expand(M) |>
    mutate( 
      se = map(.results, "se" )
    ) |>
    select( .universe, se ) |>
    unnest(se) |>
    filter(.universe == which_best) |>
    pull(se) 
  print(se_best)
  elpd_best <- multiverse::expand(M) |>
    mutate( 
      elpd = map(.results, "elpd" )
    ) |>
    select( .universe, elpd ) |>
    unnest(elpd) |>
    filter(.universe == which_best) |>
    pull(elpd)
  
  # compute the multiple comparison adjusted elpd differences
  adj_df <- multiverse::expand(M) |>
    mutate( 
      elpd = map(.results, "elpd" ),
      se = map(.results, "se" )
    ) |>
    select( .universe, elpd, se ) |>
    unnest(elpd, se) |>
    mutate(elpd_loo_diff = elpd - elpd_best,
           elpd_loo_diff_se = se_best * sqrt(1 - Rho[which_best, .universe]),
           .method = "adjusted") |>
    select(-c(elpd, se))
  
  # bind dataframes
  results <- rbind(diff_df, adj_df) |>
    filter(elpd_loo_diff_se > 0)
  return(list(df = results, which_best = which_best))
}

# compute the summary statistics from the 
results <- compute_elpd_diff(M_eeg)
diff_df <- results$df
diff_df

# fit the meta-model
hier_exec <- cmdstan_model("stan/eight_schools_rhs.stan")
stan_df <- diff_df |>
  filter(.method == "raw")
hier_data <- list(
  J = nrow(stan_df),
  y = stan_df$elpd_loo_diff,
  sigma = stan_df$elpd_loo_diff_se,
  hs_df = 1,
  hs_df_global = 1,
  hs_df_slab = 4,
  hs_scale_slab = 5,
  hs_scale_global = 2
)
hier_fit <- hier_exec$sample(data = hier_data,
                             chains = 4,
                             parallel_chains = 4,
                             refresh = 0)

# append hierarchical posterior to the results
hier_df <- as_draws_df(hier_fit$draws()) |>
  select(starts_with("theta")) |>
  summarise(
    across(
      everything(), 
      list(
        mean = ~mean(.x, na.rm = TRUE),
        se = ~sd(.x, na.rm = TRUE)
      )
    )
  ) |>
  reshape2::melt() |>
  mutate(
    model = str_split(variable, "_", simplify = T)[, 1],
    stat = str_split(variable, "_", simplify = T)[, 2]
  ) |> 
  select(!variable) |>
  reshape2::dcast(model ~ stat, value = 'value') |>
  mutate(model = as.numeric(gsub("\\D", "", model)),
         model = case_when(model >= results$which_best ~ model + 1, 
                           model < results$which_best ~ model),
         .method = "hierarchical") |>
  rename(elpd_loo_diff = mean, elpd_loo_diff_se = se, .universe = model)

## pointrange plots
## ----------------

# concatenate data into plotting data
plot_df <- rbind(diff_df, hier_df) |>
  add_row(elpd_loo_diff = 0, elpd_loo_diff_se = 0, 
          .universe = results$which_best, 
          .method = "raw") |>
  add_row(elpd_loo_diff = 0, elpd_loo_diff_se = 0, 
          .universe = results$which_best, 
          .method = "adjusted") |>
  add_row(elpd_loo_diff = 0, elpd_loo_diff_se = 0, 
          .universe = results$which_best, 
          .method = "hierarchical")

# make plot
p <- plot_df |>
  ggplot(aes(.universe, elpd_loo_diff)) +
  geom_hline(yintercept = 0, colour = "red", linetype = "dashed") +
  geom_pointrange(aes(ymin = elpd_loo_diff - elpd_loo_diff_se, 
                      ymax = elpd_loo_diff + elpd_loo_diff_se), 
                  shape = 1) +
  facet_wrap(~factor(.method, c("raw", "adjusted", "hierarchical"))) +
  xlab("Model") +
  ylab("elpdPlainHat") +
  theme_bw() + 
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        strip.background = element_blank(),
        panel.background = element_blank())
p
save_tikz_plot(p, width = 5, filename = "./tex/eeg-hierarchical.tex")




### ----

# extract the LOO-CV elpd results
( loo_df <- multiverse::expand(M_eeg) |>
    mutate( 
      elpd = map(.results, "elpd" ),
      se = map(.results, "se" )
    ) |>
    select( .universe, elpd, se ) |>
    unnest(elpd, se) )

loo_df <- multiverse::expand(M_eeg) |>
  mutate( 
    mod_loo = map(.results, "mod_loo" )
  ) |>
  select( .universe, mod_loo ) |>
  mutate(
    elpd = loo_compare(mod_loo)[,"elpd_diff"],
    se = loo_compare(mod_loo)[,"se_diff"],
  ) |>
  filter(se > 0)

# fit the meta-model
hier_exec <- cmdstan_model("stan/eight_schools.stan")
hier_data <- list(
  J = nrow(loo_df),
  y = loo_df$elpd,
  sigma = loo_df$se,
  hs_df = 1,
  hs_df_global = 1,
  hs_df_slab = 4,
  hs_scale_slab = 2,
  hs_scale_global = 1
)
hier_fit <- hier_exec$sample(data = hier_data,
                             chains = 4,
                             parallel_chains = 4,
                             refresh = 0)

# compute posterior sample statistics
( hier_diff_df <- as_draws_df(hier_fit$draws()) %>%
    select(starts_with("theta")) %>%
    summarise(
      across(
        everything(), 
        list(
          mean = ~mean(.x, na.rm = TRUE),
          sd = ~sd(.x, na.rm = TRUE)
        )
      )
    ) %>%
    reshape2::melt() %>%
    mutate(
      model = str_split(variable, "_", simplify = T)[, 1],
      stat = str_split(variable, "_", simplify = T)[, 2]
    ) %>% 
    select(!variable) %>%
    reshape2::dcast(model ~ stat, value = 'value') %>%
    mutate(model = as.numeric(gsub("\\D", "", model)), dist = "norm") %>%
    arrange(model) %>%
    mutate(model = paste0("x", as.character(model))) )

# plot Gaussian approximation to posterior
( gg_diff_dist_hier <- hier_diff_df %>% ggplot(
  aes(
    xdist = dist_normal(mean, sd),
    colour = model, 
  )
) +
    stat_slab(fill = NA) +
    scale_colour_brewer(palette = "Set2") +
    theme_classic() +
    ylab(NULL) +
    xlab(NULL) +
    ggtitle("Hierarchical model posterior") +
    theme(axis.line.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "none") )

# plot Gaussian approximation
( gg_diff_dist_analytic <- loo_df %>%
    mutate(
      dist = "norm",
      args = map2(elpd, se, list),
      model = as.factor(.universe)
    ) %>%
    ggplot(
      aes(
        xdist = dist,
        args = args, 
        colour = model
      )
    ) +
    stat_slab(fill = NA) +
    scale_colour_brewer(palette = "Set2") +
    theme_classic() +
    ylab(NULL) +
    xlab(NULL) +
    #ggtitle("LOO-CV ELPD difference") +
    ggtitle("LOO-CV elpd") +
    theme(axis.line.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "none") )

# combine plots in patchwork
x_order <- ggplot_build(gg_diff_dist_analytic)$layout$panel_scales_x[[1]]$range$range
(gg_diff_dist_analytic | gg_diff_dist_hier)
( gg_elpd_dists <- (gg_diff_dist_analytic | gg_diff_dist_hier) & 
    scale_x_continuous(limits = x_order) )

