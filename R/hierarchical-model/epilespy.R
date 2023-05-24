# load packages 
if(!requireNamespace("pacman"))install.packages("pacman")
pacman::p_load(purrr, brms, Matrix, dplyr, RColorBrewer, patchwork,
               loo, tidyr, multiverse, ggdist, ggplot2, stringr, 
               distributional, posterior, cmdstanr)

# set # of cores 
nc <- parallel::detectCores() - 1

# load data 
data("epilepsy")

# create multiverse
M_epi = multiverse()

inside(M_epi, {
  # 1. different priors: brms default and horseshoe
  # 2. include different variables: treat, zBase (interaction thereof),
  #    patient-level and visit-level random effects
  # 3. Poisson vs. negative binomial observation family 
  mod_epi <- brm(count ~
                   branch(formula, 
                          "eq_1" ~ Trt * zAge + (1 | patient),
                          #"eq_2" ~ Trt + zBase + (1 | patient),
                          #"eq_3" ~ Trt * zBase + (1 | visit),
                          #"eq_4" ~ Trt + zBase + (1 | visit),
                          "eq_5" ~ Trt * zAge + (1 | visit) + (1 | patient),
                          #"eq_6" ~ Trt + zBase + (1 | visit) + (1 | patient)
                          "eq_7" ~ Trt * zAge + (1 | visit) + (1 | patient) + (1 | obs)
                   ),
                 data = epilepsy, 
                 family = 
                   branch(family, 
                          "poisson" ~ poisson(), 
                          "negbin" ~ negbinomial()
                   ),
                 prior = 
                   branch(prior,
                          "default" ~ NULL,
                          "horseshoe" = set_prior("horseshoe(5)")
                   ),
                 cores = nc,
                 backend = 'cmdstanr',
                 silent = 2,
                 refresh = 0)
})

inside(M_epi, {
  # compute the model LOO-CV elpd
  mod_loo <- loo(mod_epi)
  elpd <- mod_loo$estimates["elpd_loo", "Estimate"]
  se <- mod_loo$estimates["elpd_loo", "SE"]
  
  # extract the pointwise results
  pointwise_loo <- mod_loo$pointwise[,"elpd_loo"]
})

# execute the multiverse analysis
execute_multiverse(M_epi)

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
    arrange(elpd, .universe) |>
    select(.universe) |>
    as.numeric()
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
  return(list(df = results, which_best = which_best, Rho = Rho))
}

# compute the summary statistics from the 
results <- compute_elpd_diff(M_epi)
diff_df <- results$df

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
  hs_scale_slab = 2,
  hs_scale_global = 1
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
plot_df |>
  ggplot(aes(.universe, elpd_loo_diff)) +
  geom_hline(yintercept = 0, colour = "red", linetype = "dashed") +
  geom_pointrange(aes(ymin = elpd_loo_diff - elpd_loo_diff_se, 
                      ymax = elpd_loo_diff + elpd_loo_diff_se), 
                  shape = 1) +
  facet_wrap(~factor(.method, c("raw", "adjusted", "hierarchical"))) +
  theme_bw() + 
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        strip.background = element_blank(),
        panel.background = element_blank())



