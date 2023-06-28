run_projpred_varsel = function(ref_model, test_data, steps = NA, ...) {
  p = ncol(test_data)-1
  n = nrow(test_data)
  if (is.na(steps)) {steps = p}
  # Run forward search
  varsel_ref = varsel(
    ref_model, 
    method='forward', 
    nterms_max=steps
  )
  
  # Do it again, but for independent test data
  d_test = list(
    data=test_data,
    offset=rep(0, n),
    weights=rep(1, n),  
    y=test_data[, 'y']
  )
  varsel_ref_test = varsel(
    ref_model, 
    method='forward', 
    nterms_max=steps, 
    d_test=d_test
  )
  
  # Combine results into df
  out = combine_varsels(
    varsel_train=varsel_ref, 
    varsel_test=varsel_ref_test
  )
  
  # Add additional information to results
  kwargs = list(...)
  for (col in names(kwargs)) {
    out[,col] = kwargs[[col]]
  }
  
  # Label base model as base
  out %>% 
    replace_na(list(solution_terms = 'base')) %>%
    tibble()
}

combine_varsels = function(varsel_train, varsel_test) {
  # Extract summary df of the varsel objects
  summary_train = summary(varsel_train, 
                          stats=c('elpd'), 
                          type = c('mean', 'se', 'diff', 'diff.se')
  )$selection
  summary_test = summary(varsel_test,
                         stats=c('elpd'), 
                         type = c('mean', 'se', 'diff', 'diff.se')
  )$selection
  # Rename columns
  names(summary_train) = c('size', 'solution_terms', 'elpd_loo', 'elpd_loo_se', 'elpd_loo_diff_ref', 'elpd_loo_diff_ref_se')
  names(summary_test) = c('size', 'solution_terms', 'elpd_test', 'elpd_test_se', 'elpd_test_diff_ref', 'elpd_test_diff_ref_se')
  
  # Merge and drop se columns
  res = merge(summary_train, summary_test, by = c('size', 'solution_terms'))
  res = res[order(res$size), ]
  rownames(res) = NULL
  
  # Add suggested size
  res$suggested_size = varsel_train$suggested_size
  
  res
}