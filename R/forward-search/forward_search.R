get_elpd_loo = function(model) {
  if(is.null(model$loo)) {
    model$loo = loo(model)
  }
  return(model$loo$estimates['elpd_loo', 'Estimate'])
}

get_elpd_test = function(model, data) {
  return(elpd(log_lik(model, newdata=data))$estimates['elpd', 'Estimate'])
}

run_forward_selection = function(model, train_data, test_data, prior, ...) {
  out = forward_selection(model, train_data, prior)
  candidates = out$candidates
  models = out$models
  # Create tibble from results
  results = map(models, function(x) {
    elpd_test_set=elpd(log_lik(x$model, newdata=test_data))
    tibble(model=list(x$model), 
           solution_terms=x$variable,
           elpd_test=elpd_test_set$estimates['elpd', 'Estimate'],
           elpd_test_se=elpd_test_set$estimates['elpd', 'SE'])
  }
  ) %>% 
    reduce(rbind) %>%
    mutate(size=row_number()-1) %>%
    rowwise() %>%
    mutate(elpd_loo=list(model$loo)) %>%
    select(!model) %>%
    ungroup()
  
  # Add additional information to results
  kwargs = list(...)
  for (col in names(kwargs)) {
    results[,col] = kwargs[[col]]
  }
  
  # This adds the loo differences to the output tibble for convinience
  results = compute_elpd_diff(results)
  return(list(results = results, candidates = candidates))
}

forward_selection = function(model, train_data, prior) {
  # Initial set of predictors
  cols = names(train_data)
  predictors = cols[which('y'!=cols)]
  steps = length(predictors)
  # Forward step function returns this format
  # Field `variable` contains the choses predictor
  new = list(model=model, variable=model$variable)
  
  # Collect forward selection path here
  models = list()
  candidates = list()
  
  models[[1]] = new
  k = 1
  while(k <= steps) {
    # Choose the next variable
    out = forward_step(new$model, train_data, predictors, prior)
    new = out$selected
    # Drop the selected variable from candidates
    idx = which(new$variable != predictors)
    predictors = predictors[idx]
    models[[k+1]] = new
    candidates[[k]] = out$candidates
    k = k + 1
  }
  return(list(models=models, candidates=candidates))
}


update_model = function(model, train_data, x, prior) {
  formula_new = update(model$formula, paste0('~ . + ', x))
  n = nrow(train_data)
  if (class(model) == 'brmsfit') {
    model_new = update(model, 
                       newdata=train_data, 
                       formula=formula_new,
                       prior=prior,
                       refresh=0)
  } else {
    stop()
  }
  model_new$loo = loo(model_new)
  return(list(model=model_new, variable=x))
}

forward_step = function(model, train_data, predictors, prior) {
  candidate_models = lapply(predictors, function(x) update_model(model, train_data, x, prior))
  selected = NULL
  current_loo = -10e6
  for (candidate in candidate_models) {
    new_loo = get_elpd_loo(candidate$model)
    if (new_loo > current_loo) {
      current_loo = new_loo
      selected = candidate
    }
  }
  # Extract ELPD-LOOs of the candidate models
  candidate_elpd_loo = candidate_models %>% 
    map('model') %>% 
    map('loo') %>% 
    map_dbl(function(x) x$estimates['elpd_loo', 'Estimate'])
  candidate_elpd_loo_se = candidate_models %>% 
    map('model') %>% 
    map('loo') %>% 
    map_dbl(function(x) x$estimates['elpd_loo', 'SE'])
  candidate_var = candidate_models %>% map_chr('variable')
  
  candidates = candidate_models %>% 
    map('model') %>% 
    map(function(x) loo_diff(x$loo, model$loo)) %>%
    map(as_tibble) %>%
    reduce(rbind) %>% 
    mutate(candidate_elpd_loo, candidate_elpd_loo_se, candidate_var)
    
  return(list(selected=selected, candidates=candidates))
}


compute_elpd_diff = function(results) {
  elpd_diffs = map2(results$elpd_loo, lag(results$elpd_loo), loo_diff)
  elpd_diffs = do.call(rbind, elpd_diffs) %>% 
    as_tibble() %>% 
    unnest(cols = c(elpd_loo_diff, elpd_loo_diff_se))
  results$elpd_loo_diff = elpd_diffs$elpd_loo_diff
  results$elpd_loo_diff_se = elpd_diffs$elpd_loo_diff_se
  results$elpd_loo_se = results$elpd_loo %>% 
    map_dbl(function(x) x$estimate['elpd_loo', 'SE'])
  results$elpd_loo = results$elpd_loo %>% 
    map_dbl(function(x) x$estimate['elpd_loo', 'Estimate'])
  results
}


loo_diff = function(loo1, loo2) {
  if(is.null(loo1)) return(NA)
  if(is.null(loo2)) return(NA)
  stopifnot(nrow(loo1$pointwise) == nrow(loo2$pointwise))
  
  n = nrow(loo1$pointwise)
  diff_se = sd(loo1$pointwise[, 'elpd_loo'] - loo2$pointwise[, 'elpd_loo']) * sqrt(n)
  diff = sum(loo1$pointwise[, 'elpd_loo'] - loo2$pointwise[, 'elpd_loo'])
  return(list(elpd_loo_diff = diff, elpd_loo_diff_se = diff_se))
}