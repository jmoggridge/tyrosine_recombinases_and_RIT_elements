library(tidyverse)
library(tidymodels)
library(workflowsets)

set.seed(1)

# make a small dataset with stratified sampling for code development
df <- 
  read_rds('./data/smart_refseqs_hmm_scores.rds') |>
  select(-id) |> 
  relocate(subfamily, contains('prot_')) |> 
  group_by(subfamily) |> 
  sample_frac(0.1) |> 
  ungroup() |> 
  mutate(across(Arch1:Xer, ~replace_na(.x, 0))) |> 
  janitor::clean_names()

# count kmers


glimpse(df)
summary(df)
summary(df$subfamily)


library(corrr)

cor_df <- df |> 
  select(arch1:xer) |> 
  mutate(across(everything(), ~(.x - mean(.x)) / sd(.x))) |> 
  correlate() |> 
  pivot_longer(-term, names_to = 'term2', values_to = 'cor') |> 
  mutate(cor = ifelse(is.na(cor), 0, cor)) |> 
  arrange(-cor) |> 
  print(n=50)

cor_df |> 
  group_by(term) |> 
  summarize(mean = mean(cor)) |> 
  arrange(-mean)


# set up models
model_ls <- list(
  
  # random forest
  rf_mod <- 
    rand_forest(mtry = tune(), min_n = tune(), trees = tune()) |> 
    set_engine('ranger') |> 
    set_mode('classification'),
  
  # k-nearest neighbors
  knn_mod <- 
    nearest_neighbor(
      neighbors = tune(), 
      weight_func = tune(), 
      dist_power = tune()) |> 
    set_engine('kknn') |> 
    set_mode('classification'),
  
  # xgboost
  xgb_mod <- boost_tree(mtry = tune(), trees = tune(), min_n = tune()) |> 
    set_mode('classification') |> 
    set_engine('xgboost')
  # logistic regression
  # neural network
  # 
)
model_ls

# set up performance measures
metrix <- metric_set(precision, recall, f_meas, kap, bal_accuracy, mcc)

# set up resampling scheme
df_resampling <- 
  nested_cv(
    df, 
    outside = vfold_cv(v = 5, strata = subfamily),
    inside = bootstraps(times = 25, strata = subfamily)
    )

# data preprocessing
df_recipe <- recipe(subfamily ~ ., data = df) |> 
  step_normalize(all_numeric_predictors()) |> 
  update_role(c(contains('prot_')), new_role = "sample id") 

df_recipe |> prep()
df_recipe |> prep() |>  juice()

# cross recipe with models to create multiple workflows
model_wf <- workflow_set(preproc = list(df = df_recipe), models = model_ls)

# tuning parameter search with bayesian optimization
tune_rs <- 
  model_wf |> 
  workflow_map(
    fn = 'tune_bayes', 
    resamples = df_resampling,
    metrics = metrix,
    verbose = T)
# 
# 
# tune_bayes(
#   object,
#   preprocessor,
#   resamples,
#   ...,
#   iter = 10,
#   param_info = NULL,
#   metrics = NULL,
#   objective = exp_improve(),
#   initial = 5,
#   control = control_bayes()
# )




