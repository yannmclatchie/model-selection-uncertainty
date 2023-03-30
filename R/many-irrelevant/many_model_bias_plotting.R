library(ggplot2)
library(geomtextpath)
library(dplyr)

# read in results
results <- vroom("") # TODO

# plot results
plotting_df <- results %>% group_by(beta, iter) %>%
  mutate(best_loo = max(loo_diff)) %>%
  group_by(beta) %>%
  mutate(mean_best_loo = mean(best_loo)) %>%
  right_join(out %>% 
               filter(model==1) %>%
               group_by(beta, iter) %>%
               mutate(true_diff = mean(test_diff)) %>%
               group_by(beta) %>%
               mutate(mean_true_diff = mean(test_diff))) %>%
  select(c("beta", "K", "mean_best_loo", "mean_true_diff")) %>%
  unique() %>%
  reshape2::melt(id = c("beta", "K"))
( p <- plotting_df %>% 
    ggplot(aes(x = beta, y = value, colour = variable, label = variable)) + 
    geom_hline(yintercept=0, colour = "grey") + 
    geom_vline(xintercept=0, colour = "grey") + 
    geom_point() +
    geom_labelsmooth(aes(hjust = variable), method = glm, 
                     se = FALSE, linetype = 3) +
    facet_wrap(vars(K)) +
    ylab("elpd") + 
    xlab("beta") + 
    theme_classic() +
    theme(legend.position="none") )