library(tidyverse)
library(tidymodels)
library(workflowsets)
library(themis)
library(doParallel)
library(doFuture)

full_data <- read_rds('./data/full_classifier_dataset.rds')

# downsample to max 500 rows per class
set.seed(1)
df <- full_data |> 
  select(-c(prot_seq, prot_description)) |> 
  group_by(subfamily) |> 
  slice_sample(n=500) |> 
  ungroup() |> 
  mutate(across(Arch1:Xer, ~replace_na(.x, 0)))
summary(df$subfamily)

rm(full_data)


# split data 75 % training / 25 % final test
init_split <- initial_split(df, strata = subfamily, pool = 0.02)
train_df <- training(init_split)
test_df <- testing(init_split)

train_df |> count(subfamily) |> print.data.frame()

# data preprocessing spec to perform within cv
recip <- recipe(subfamily ~ ., data = df) |> 
  update_role(c(contains('prot_')), new_role = "sample id") |> 
  step_zv(all_numeric_predictors()) |> 
  step_normalize(all_numeric_predictors()) |> 
  step_smote(all_outcomes())
recip

recip |> prep() |> juice() |> count(subfamily) |> print.data.frame()

### models
# svm_mod <- svm_rbf(mode = "classification", cost = tune()) %>% 
#   set_engine("kernlab")
# svm_mod

rf_mod <- rand_forest(mtry = tune(), min_n = tune()) |> 
  set_mode('classification') |> 
  set_engine('ranger')
rf_mod

# combine recipe and model into workflow
wf <- workflow() |> add_recipe(recip) |> add_model(rf_mod)
wf
tune_args(wf)

rf_set <- 
  parameters(wf) %>% 
  update(mtry = mtry(range = c(1L, 420L)))
rf_set

# set up performance measures
# perf_metrics <- metric_set(mcc, precision, recall, f_meas, kap, bal_accuracy, )
# perf_metrics

# resampling procedure for training data
results <- nested_cv(
  train_df, 
  outside = vfold_cv(strata = 'subfamily', repeats = 3, v = 3, pool = 0.01),
  inside = vfold_cv(strata = 'subfamily', repeats = 3, v = 3, pool = 0.01) 
  )
results

# repeated vfold cv
# resamp <- vfold_cv(train_df, repeats = 3, v = 3, strata = subfamily, pool = 0.01)
# 
# class(results$inner_resamples[[1]])
# results$inner_resamples[[1]]
# results$inner_resamples[[1]]$splits[[1]]
# analysis(results$inner_resamples[[1]]$splits[[1]]) |> 
#   count(subfamily) |> print.data.frame()
# assessment(results$inner_resamples[[1]]$splits[[1]]) |> 
#   count(subfamily) |> print.data.frame()
# summary(results$inner_resamples[[2]]$splits[[1]]$data$subfamily)

# setup parallelization
all_cores <- parallel::detectCores(logical = FALSE)
registerDoFuture()
cl <- makeCluster(all_cores)
plan(cluster, workers = cl)

rs <- tune_bayes(
  wf, 
  resamples = results$inner_resamples[[1]],
  initial = 5,
  iter = 20,
  param_info = rf_set,
  # How to measure performance?
  # metrics = metric_set(mcc),
  control = control_bayes(no_improve = 10, verbose = TRUE)
  )

rs |> collect_metrics()
rsmet <- rs |> 
  show_best(metric = 'roc_auc') |> print.data.frame()
rsmet


# map tune_bayes over inner resamples??

write_rds(rs, './data/tunebayes_rs.rds')


# `object` will be an `rsplit` object from our `results` tibble
# `cost` is the tuning parameter

# evaluate_model_rmse <- function(object, cost = 1) {
#   y_col <- ncol(object$data)
#   mod <- 
#     svm_rbf(mode = "regression", cost = cost) %>% 
#     set_engine("kernlab") %>% 
#     fit(y ~ ., data = analysis(object))
#   
#   holdout_pred <- 
#     predict(mod, assessment(object) %>% dplyr::select(-y)) %>% 
#     bind_cols(assessment(object) %>% dplyr::select(y))
#   rmse(holdout_pred, truth = y, estimate = .pred)$.estimate
# }





# # models
# 
# # set up models
# model_ls <- list(
#   # random forest
#   rf_mod <- 
#     rand_forest(mtry = tune(), min_n = tune(), trees = tune()) |> 
#     set_engine('ranger') |> 
#     set_mode('classification'),
#   # k-nearest neighbors
#   knn_mod <- 
#     nearest_neighbor(neighbors = tune(), weight_func = tune(), dist_power = tune()) |> 
#     set_engine('kknn') |> 
#     set_mode('classification'),
#   # svm non-linear kernel
#   svm_mod <- svm_rbf(mode = "classification", cost = tune()) %>% 
#     set_engine("kernlab")
#   
#   # xgboost
#   # xgb_mod <- boost_tree(mtry = tune(), trees = tune(), min_n = tune()) |> 
#   #   set_mode('classification') |> 
#   #   set_engine('xgboost')
#   # logistic_reg
#   
#   # neural network
#   # ...
# )\



# # tuning parameter search with bayesian optimization
# tune_rs <- 
#   model_wf |> 
#   workflow_map(
#     fn = 'tune_bayes', 
#     resamples = df_resampling,
#     metrics = metrix,
#     verbose = T)

# df_recipe |> prep()
# df_recipe |> prep() |>  juice()
# 
# # cross recipe with models to create multiple workflows
# model_wf <- workflow_set(preproc = list(df = df_recipe), models = model_ls)
# 
