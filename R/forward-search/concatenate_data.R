library(tidyverse)

read_data = function(what, folder) {
  stopifnot(what %in% c('projpred', 'forward_search', 'forward_search_normal', 'candidates'))
  path = file.path('results/forward_variable_selection', folder)
  files = dir(path)
  res = list()
  for (file in files) {
    d = readRDS(file.path(path, file))
    
    if(what == 'candidates') {
      candidates = d[[what]]
      fs = d[['forward_search']]
      n = unique(fs$n)[1]
      rho = unique(fs$rho)[1]
      iter = unique(fs$iter)[1]
      candidates = map2(candidates, seq_along(candidates)-1, function(x, y) {
        x$size = y
        x$n = n
        x$rho = rho
        x$iter = iter
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

fs <- read_data('forward_search', 'data/results/forward-search') |>
  write_csv('data/results/forward-search/forward_search.csv')

candidates_fs = read_data('candidates', 'data/results/forward-search') |>
  inner_join(fs %>% select(n, rho, iter, size, elpd_loo, elpd_loo_se)) |>
  write_csv('data/results/forward-search/candidates.csv')


