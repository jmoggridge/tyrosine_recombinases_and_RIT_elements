# name for directory in project folder to store nested CV files
run_name <- 'test_debug'

source('./R/00_functions.R')

# create directory structure for classifier files (alignments, hmms, results for each resample)
out_path <- glue(here::here(), '/', run_name)
system(glue('mkdir {out_path}'))

list('align','hmm', 'hmmsearch', 'results') |> 
  map(~glue(out_path, '/', .x)) |>
  map(~system(glue('mkdir ', .x)))

out_path

train <- read_rds('./data/classif_train_set.rds')
train |> ggplot(aes(subfamily)) + geom_bar()


## test prep functions ----

# fold <- train |> 
#   group_by(subfamily) |> 
#   slice_sample(n = 500, replace = F) |> 
#   ungroup()
# fold |> ggplot(aes(subfamily)) + geom_bar()
# 
# sum(is.na(fold))
# fold |> filter(is.na(dom_seq)) |> count(subfamily)
# prep_domains_df(fold)
# 
# split_id = 'test1'
# cat('\n', green$bold(glue("Preparing fold: {split_id}")))
# # create directories
# system(glue('mkdir {out_path}/align/{split_id}'))
# system(glue('mkdir {out_path}/hmm/{split_id}'))
# system(glue('mkdir {out_path}/hmmsearch/{split_id}'))
# 
# # check stratification is adequate / no empty classes
# strat_test <- all(
#   fold |> count(subfamily) |> pull(n) %>% all(.data > 10), 
#   fold |> count(subfamily) |> pull(n) %>% all(.data > 10)
# )
# cat('\n', white(glue("Train and test split contain > 10 obs per subfamily? {strat_test}")))
# 
# 
# 
# 
# prep <-  fold |> prep_domains_df()
# prep |> ggplot(aes(fct_rev(subfamily))) + geom_bar() + coord_flip()
# 
# als <- build_alignments_library(domains = prep, split_id = split_id, out_path = out_path)
# als
# 
# 
# # train hmm
# cat(white(' Building HMMs...'))
# hmm_check <- build_hmm_library(split_id = split_id, out_path = out_path)
# hmmbuild_test <- all(hmm_check$hmmbuild == 0)
# cat('\n', white(glue('hmmbuild calls all finished without error? {hmmbuild_test}')))
# 
# # hmm scores test and train
# cat('\n', white('Scoring sequences against HMM library...'))
# tic()
# train_search <- 
#   fold |> 
#   hmmsearch_scores(split_id = split_id, tag = 'train', out_path = out_path)
# 
# hmmsearch_check <- all(train_search$hmmsearches == 0)
# cat(white(glue(' HMM searches all finished without error? {hmmsearch_check}')))
# 
# # parse hmmscores and add to data frame
# train <-  
#   fold |> 
#   join_hmmsearches(files = train_search$out_path)
# train



## ------



# `NA` values are not allowed when using `step_nearmiss`

resamples <- nest_cv$inner_resamples[[1]]
resamples


# do_inner_cv <- function(resamples, out_path){

cat("\n", blue$bold("Starting inner CV"), '\n')
# prepare scored datasets for each fold
prep <- resamples |> 
  mutate(prep = map2(inner_splits, inner_id, 
                     ~ prep_data(.x, .y, out_path = out_path))) |>
  unnest(prep)

prep |> select(train) |> unnest(c(train)) |> is.na() |> sum()

map(names(prep |>  select(train) |> unnest(c(train)) ),
    ~{prep |> select(train) |> unnest() |> filter(is.na({{.x}}))})

beepr::beep()

models
# 
# model <- models |>  pull(spec) |> pluck(120)
# model
# train <- bind_rows(prep$train) |> 
#   select(-contains('seq'))
# test <- prep$test[[1]] |> 
#   select(-contains('seq'))
# 
# train
# test
# # |>
#   # ungroup()
# sum(is.na(train))
# sum(is.na(test))

# 
# ## Fit + evaluate model on single fold
#   # eval_model <- function(model, train, test){
#     # specify evaluation functions
# class_metrics <-
#   metric_set(mcc, kap, sensitivity, specificity, precision, recall,
#              f_meas, ppv, npv, accuracy, bal_accuracy)
# 
# # specify formula and scaling recipe
# recip <-
#   recipe(subfamily ~ ., data = train) |>
#   update_role(acc, new_role = "id") |>
#   step_downsample(subfamily, under_ratio = 10) |>
#   step_smote(subfamily, over_ratio = 1) |>
#   step_normalize(all_predictors())
# 
# 
# # create model workflow and fit to training data
# fit_wf <-
#   workflow() |>
#   add_model(model) |>
#   add_recipe(recip) |>
#   fit(data = train)
# 
# # predict hold out set
# predictions <- fit_wf |>
#   predict(test |> select(-subfamily)) %>%
#   bind_cols(test)
# # collect metrics
# predictions |>
#   class_metrics(truth = subfamily, estimate = .pred_class)
#   # }
# 
# predictions |> conf_mat(truth = subfamily, estimate = .pred_class)

# 
train <- prep$train[[1]]
test <- prep$test[[1]]
# rs1 = map2(train1, test1, ~ eval_model_set(models, .x, .y))

model = models$spec[[3]]

## Fit + evaluate model on single fold
# eval_model <- function(model, train, test){
  
  train <- train |> select(-contains('seq'))
  test <- test |> select(-contains('seq'))
  # specify evaluation functions
  class_metrics <- 
    metric_set(mcc, kap, sensitivity, specificity, precision, recall, 
               f_meas, ppv, npv, accuracy, bal_accuracy)
  
  # specify formula and scaling recipe
  recip <- 
    recipe(subfamily ~ ., data = train) |> 
    update_role(acc, new_role = "id") |> 
    step_downsample(subfamily, under_ratio = 3) |> 
    step_smote(subfamily, over_ratio = 0.5) |> 
    step_normalize(all_predictors())
  
  # create model workflow and fit to training data
  fit_wf <- 
    workflow() |>
    add_model(model) |> 
    add_recipe(recip) |>
    fit(data = train)
  
  # predict hold out set
  predictions <- fit_wf |> 
    predict(test |> select(-subfamily)) %>% 
    bind_cols(test)
  # collect metrics
  predictions |> 
    class_metrics(truth = subfamily, estimate = .pred_class)
}

# 
# train2 <- prep$train[[2]]
# test2 <- prep$test[[2]]
# eval_model_set(models = models, train = train2, test = test2) |> View()
# 
train3 <- prep$train[[3]]
test3 <- prep$test[[3]]
eval_model_set(models = models, train = train3, test = test3) |> View()

train <- prep$train[[4]]
test <- prep$test[[4]]
eval_model_set(models = models, train = train4, test = test4) |> View()


# 
# cat('\n', blue$bold("Evaluating models"))
# tic()
# 
# models <- models[1:5,]
# 
# eval_model_set(models = models, train = train2, test = test2) |> View()

# evaluate model set on each fold
fold_res <- prep |>
  mutate(rs = map2(train, test,
                   ~ eval_model_set(models = models, train = .x, test = .y)))

# |>
#   unnest(rs)
# toc()
# 
# 

