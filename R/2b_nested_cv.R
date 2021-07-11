## Classifier training / assessment by repeated nested-stratified CV

# Performs a nested CV where the training data are used to create domain subfamily alignments and train HMMs. Then the HMMs aer used to score train and test sequences. The scores are used to predict the class of the test data, and performance metrics for the different models are compared in the inner CV (model tuning) and outer CV (model selection).

## done: TODO selection of models
## done: TODO threshold classification
## done: TODO fix splits labels
## done: TODO selection of k and rep for outer and inner CV: 3x3

## Libraries -------------------------------------------------------------

library(tidyverse)
library(tidymodels)
library(themis)
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

# unfitted model specifications
models <- read_rds('./data/unfitted_parsnip_model_set.rds')

set.seed(123)

## Training data
# do data splitting first './R/2a_data_splitting.R'
train <- read_rds('./data/classif_train_set.rds')

# train |> ggplot(aes(fct_rev(subfamily))) + geom_bar() + coord_flip()

# # TODO remove downsampling
# train <- train |> 
#   group_by(subfamily) |> 
#   slice_sample(n = 300, replace = F) |> 
#   ungroup()
#   
#   
  
  
## Directories --------------------------------------------------------

# name for directory in project folder to store nested CV files
run_name <- '3x3-fold_07-08'


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

nest_cv |> unnest(inner_resamples) |> pull(inner_splits) |> pluck(1) |> analysis() |> 
  count(subfamily) |> pull(n) |> maxmin(
  )
maxmin <- function(x)  max(x)/min(x)

rm(train)

# MAIN / TRAIN & TEST ----------------------------------------

# **evaluate models with nested CV** all functions are in 00_functions.R
tic()
nest_cv_results <-  fit_nested_cv(
  nestcv = nest_cv, 
  models = models,
  out_path = out_path)
toc()
beepr::beep()

write_rds(nest_cv_results, 
          glue('./results/{run_name}_nest_cv_results.rds'), 
          compress = 'gz')

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

nest_cv_summary |> print(n=100)

# evaluate maxscore/threshold classifications for outer folds and summarize
thresh_res <- nest_cv_results |> 
  mutate(thresh_res = map2(train, test,
                           ~eval_threshold_classififer(.x, .y))) |> 
  select(outer_id, thresh_res) |> 
  unnest(cols = c(thresh_res)) |> 
  group_by(.metric) |> 
  summarise(model_type = 'maxscore + threshold',
            mean = mean(.estimate, na.rm = T),
            err = sd(.estimate, na.rm = T), 
            n_folds = sum(!is.na(.estimate)),
            values = list(.estimate),
            ) |> 
  relocate(model_type)
thresh_res

# combine model summary and rule-based classifier summary
final_summary <- bind_rows(nest_cv_summary, thresh_res) 

final_summary |> print(n=200)

write_rds(final_summary, glue('./results/{run_name}_nest_cv_summary.rds'))

