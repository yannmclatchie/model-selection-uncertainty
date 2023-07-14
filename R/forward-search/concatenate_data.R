library(tidyverse)

# read in files
files <- list.files("data/results/forward-search/", full.names = TRUE)

# read and write R2D2 forward search results
out_r2d2 <- files |> 
  map(readRDS) |> 
  map(rbind) |> 
  map(\(x) x[,"forward_search"]$forward_search) |> 
  bind_rows() |>
  write_csv("data/results/forward-search/fs.csv")

# read and write normal forward search results
files |> 
  map(readRDS) |> 
  map(rbind) |> 
  map(\(x) x[,"forward_search_normal"]$forward_search_normal) |> 
  bind_rows() |>
  write_csv("data/results/forward-search/fs_normal.csv")

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
      rho = unique(fs$rho)[1]
      iter = unique(fs$iter)[1]
      name = unique(fs$name)[1]
      candidates = map2(candidates, seq_along(candidates)-1, function(x, y) {
        x$size = y
        x$n = n
        x$rho = rho
        x$iter = iter
        x$name = name
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

# read and write R2D2 candidates
read_data('candidates', 'data/results/forward-search') |>
  write_csv('data/results/candidates_fs.csv')

# read and write normal candidates
read_data('candidates_normal', 'data/results/forward-search') |>
  write_csv('data/results/candidates_fs_normal.csv')
