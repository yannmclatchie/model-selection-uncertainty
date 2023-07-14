library(readr)

train_df <- read_csv("data/datasets/SPECTF/spect_train.csv")
test_df <- read_csv("data/datasets/SPECTF/spect_test.csv")
df <- rbind(train_df, test_df) |> rename(y = OVERALL_DIAGNOSIS)
df
sum(is.na(df))


prior = c(
  prior(horseshoe(df = 3, scale_slab = 50, par_ratio = 0.3), class = 'b'),
  prior('normal(0, 2.5)', class='Intercept')
)
reference = brm(y ~ ., 
                data = df, 
                prior = prior, 
                family = bernoulli(),
                backend = 'cmdstanr',
                stan_model_args=list(stanc_options = list("O1")),
                control=list(adapt_delta=0.95),
                silent=2, 
                refresh=0)
loo(reference)

# generate data
k <- 10 # number of folds

# Generating random indices 
id <- sample(rep(seq_len(k), length.out=nrow(df)))
table(id)

# lapply over them:
indicies <- lapply(seq_len(k), function(a) list(
  n_test = sum(id==a),
  n_train = sum(id!=a),
  train_df = df[id!=a, ],
  test_df = df[id==a, ],
  data = "heart"
))
indicies

saveRDS(indicies, file = "model-selection-uncertainty/data/datasets/heart.RDS")

data <- readRDS("model-selection-uncertainty/data/datasets/heart.RDS")
cv_fold <- 2
current_data <- data[[cv_fold]]

df <- current_data$train_df
df_test <- current_data$test_df
n <- current_data$n_train
n_test <- current_data$n_test
data_name <- current_data$data

print(dim(df))
print(dim(df_test))
print(n)
print(n_test)
print(data_name)
