## Classifier training / assessment by repeated nested-stratified CV

# Performs a nested CV where the training data are used to create domain subfamily alignments and train HMMs. Then the HMMs aer used to score train and test sequences. The scores are used to predict the class of the test data, and performance metrics for the different models are compared in the inner CV (model tuning) and outer CV (model selection).

## done: TODO selection of models
## done: TODO threshold classification
## done: TODO fix splits labels
## done: TODO selection of k and rep for outer and inner CV: 3x3

## Libraries -----------------------------------------------------------------

library(tidyverse)
library(tidymodels)
library(Biostrings)
library(future)
library(furrr)
library(DECIPHER)
library(glue)
library(here)
library(tictoc)
library(crayon)
library(beepr)
# option to get `future::plan()` to work on my macbook
options(parallelly.makeNodePSOCK.setup_strategy = "sequential")
# functions 
source('./R/00_functions.R')
# models
source('./R/00_get_model_specs.R')

set.seed(123)

## Training data
# do data splitting first
# source('./R/2a_data_splitting.R')
train <- read_rds('./data/classif_train_set.rds')



## Directories ----------------------------------------------------------------

# name for directory in project folder to store nested CV files
run_name <- '3x3-fold_06-30'
# TODO transfer files to nestcv_06-30 folder

# create directory structure for classifier files (alignments, hmms, results for each resample)
out_path <- glue(here::here(), '/', run_name)
system(glue('mkdir {out_path}'))

list('align','hmm', 'hmmsearch', 'results') |> 
  map(~glue(out_path, '/', .x)) |>
  map(~system(glue('mkdir ', .x)))

## Setup Nested CV --------------------------------------------------

nest_cv <- 
  nested_cv(
    train, 
    outside = vfold_cv(v = 3, repeats = 3, strata = subfamily),
    inside = vfold_cv(v = 3, repeats = 3, strata = subfamily)
  ) |> 
  transmute(
    outer_id = glue('out_{id}_{id2}') |> str_remove_all('Repeat|Fold'),
    outer_splits = splits,
    inner_resamples
  ) |>
  unnest(inner_resamples) |> 
  transmute(
    outer_id, 
    outer_splits,
    inner_id = glue('{outer_id}_in_{id}_{id2}'),
    inner_splits = splits
  ) |> 
  nest(inner_resamples = c(inner_id, inner_splits))

nest_cv |> unnest(inner_resamples)

rm(train)

# MAIN / TRAIN & TEST ----------------------------------------

# **evaluate models with nested CV** all functions are in 00_functions.R
tic()
nest_cv_results <-  fit_nested_cv(nest_cv)
toc()
beepr::beep()

write_rds(nest_cv_results, glue('./results/{run_name}_results.rds'))

# collect results and get means
nest_cv_summary <- 
  nest_cv_results |> 
  select(results) |> 
  unnest(cols = c(results)) |> 
  group_by(model_type, .metric) |> 
  summarize(
    mean = mean(.estimate, na.rm = T),
    err = sd(.estimate, na.rm = T), 
    n_folds = sum(!is.na(.estimate)),
    values = list(.estimate),
    .groups = 'drop'
  ) 

# evaluate maxscore/threshold classifications for outer folds and summarize
thresh_res <- nest_cv_results |> 
  mutate(thresh_res = map2(train, test, ~eval_threshold_classififer(.x, .y))) |> 
  select(outer_id, thresh_res) |> 
  unnest(cols = c(thresh_res)) |> 
  group_by(.metric) |> 
  summarise(model_type = 'maxscore_threshold',
            mean = mean(.estimate, na.rm = T),
            err = sd(.estimate, na.rm = T), 
            n_folds = sum(!is.na(.estimate)),
            values = list(.estimate),
            ) |> 
  relocate(model_type)

# combine model summary and rule-based classifier summary
final_summary <- bind_rows(nest_cv_summary, thresh_res) 

write_rds(final_summary, glue('./results/{run_name}_result_summary.rds'))

final_summary |>
  filter(.metric %in% c('mcc', 'bal_accuracy', 'f_meas', 'kap',
                       'precision', 'recall','sens', 'spec')) |>
  View()
  

library(gt)

# compare models by metrics
final_summary |>
  select(-values, -n_folds) |> 
  mutate(mean = str_extract(as.character(mean), '......'),
         err = str_extract(as.character(err), '......')) |> 
  transmute(
    model = str_replace_all(model_type, '_', ' '),
    metric = str_replace_all(.metric, '_', ' '),
    mean_sd = glue('{mean} ± {err}')) |> 
  filter(metric %in% c('mcc', 'bal_accuracy', 'f_meas', 'kap',
                       'precision', 'recall','sens', 'spec')) |> 
  pivot_wider(id_cols = model, names_from = metric, values_from = mean_sd) |> 
  gt() |> 
  tab_header(title = "Performance metrics from nested 3-fold cross-validation repeated 3 times") 


# check out predictions and see which classes are misclassified most often


# sample of HMM scores
inner_cv_mcc <- 
  nest_cv_results |> 
  select(outer_id, inner_cv) |> 
  unnest(inner_cv) |>
  filter(.metric == 'mcc') |> 
  left_join(models, by = c("model_type", "model_id"))






