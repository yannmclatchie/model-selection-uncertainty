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
                          "eq_1" ~ zBase + (1 | patient) + (1 | visit),
                          #"eq_2" ~ Trt * zBase + (1 | patient) + (1 | visit),
                          "eq_3" ~ zBase + (1 | patient),
                          #"eq_4" ~ Trt * zBase + (1 | patient),
                          "eq_5" ~ Trt + (1 | patient) + (1 | visit),
                          "eq_7" ~ Trt + (1 | patient),
                          #"eq_8" ~ Trt + zBase + (1 | patient) + (1 | obs),
                          #"eq_9" ~ Trt + (1 | patient) + (1 | visit) + (1 | obs),
                          #"eq_10" ~ Trt + zAge + (1 | patient) + (1 | visit) + (1 | obs),
                          "eq_11" ~ Trt * zBase + (1 | patient) + (1 | visit) + (1 | obs)
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
})

# execute the multiverse analysis
execute_multiverse(M_epi)

# extract the LOO-CV elpd results
( loo_df <- multiverse::expand(M_epi) |>
  mutate( 
    elpd = map(.results, "elpd" ),
    se = map(.results, "se" )
  ) |>
  select( .universe, elpd, se ) |>
  unnest(elpd, se) )

loo_df <- multiverse::expand(M_epi) |>
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
hier_exec <- cmdstan_model("stan/eight_schools_rhs.stan")
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
    ggtitle("LOO-CV ELPD difference") +
    theme(axis.line.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "none") )

# combine plots in patchwork
x_order <- ggplot_build(gg_diff_dist_analytic)$layout$panel_scales_x[[1]]$range$range
(gg_diff_dist_analytic | gg_diff_dist_hier)
( gg_elpd_dists <- (gg_diff_dist_analytic | gg_diff_dist_hier) & 
    scale_x_continuous(limits = x_order) )

