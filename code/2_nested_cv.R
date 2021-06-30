## Classifier training / assessment by repeated nested-stratified CV

# Performs a nested CV where the training data are used to create domain subfamily alignments and train HMMs. Then the HMMs aer used to score train and test sequences. The scores are used to predict the class of the test data, and performance metrics for the different models are compared in the inner CV (model tuning) and outer CV (model selection).

# TODO selection of models
# TODO selection of k and rep for outer and inner CV
# TODO threshold classification
# TODO fix splits labels

## Libraries -----------------------------------------------------------------

library(tidyverse)
library(tidymodels)
library(Biostrings)
library(furrr)
library(DECIPHER)
library(glue)
library(here)
library(tictoc)
library(crayon)

# functions 
source('./code/00_functions.R')
 

## Directories ----------------------------------------------------------------

# name for directory in project folder to store nested CV files
dir <- 'classifier_01/'

# create directory structure for classifier files (alignments, hmms, results for each resample)
out_path <- glue(here::here(), '/', dir)
system(glue('mkdir {out_path}'))

list('align','hmm', 'hmmsearch', 'results') |> 
  map(~glue(here::here(), '/', dir, .x)) |>
  map(~system(glue('mkdir {.x}')))
rm(dir)


### Models  -------------------------------------------------------------------

# update model specifications using parameters grid
make_models <- function(model, grid, name){
  models <- map(
    seq_len(nrow(grid)), ~{
      mods <- model |> 
        update(grid[.x, ])
    })
  tibble(spec = models) |> 
    mutate(model_type = name, 
           model_id = row_number()) |> 
    bind_cols(grid) |> 
    nest(params = names(grid))
}

# decision tree models
tree_models <- 
  decision_tree(mode = 'classification') |> 
  set_engine('rpart') |> 
  make_models(
    name = 'decision_tree_rpart',
    grid = grid_max_entropy(
      size = 20, tree_depth(), cost_complexity(), min_n(range = c(2L, 40L))
    ))

# logistic reg models
log_reg_models <- 
  multinom_reg(mode = 'classification') |> 
  set_engine('glmnet') |> 
  make_models(name = 'multinom_reg_glmnet',
              grid = grid_max_entropy(mixture(), penalty(), size = 20))
knn_models <- 
  nearest_neighbor(mode = 'classification') |> 
  set_engine('kknn') |> 
  make_models(name = 'nearest_neighbor',
              grid = grid_regular(neighbors(), levels = 20)
  )
# combine all model specs
models <- bind_rows(tree_models, log_reg_models, knn_models)

rm(tree_models, log_reg_models, knn_models)


# MAIN --------------------------------------------------------------------

## Combine datasets  ------

# integrase dataset
smart_df <- 
  read_rds('./data/SMART/smart_df.rds') |>
  select(subfamily, acc, description, prot_seq, dom_seq) |> 
  mutate(subfamily = as_factor(subfamily))

# 20 subfamilies for classification
subfamilies <- unique(smart_df$subfamily)

# negative examples (non-integrases)
non_integrases <- 
  read_rds('./data/non_integrase_seqs/nonint_df.rds') |> 
  mutate(dom_seq = NA,
         subfamily = as_factor(subfamily)) |> 
  select(subfamily, acc, description, prot_seq, dom_seq)

set.seed(123)
# combine datasets
full_dataset <- 
  bind_rows(smart_df, non_integrases) |> 
  # TODO (remove downsampling)
  # downsample to 10000
  group_by(subfamily) |> 
  slice_sample(n = 15000, replace = F) |> 
  ungroup()

rm(smart_df, non_integrases)

## Initial split -------------------------------------------------

set.seed(54321)

df_split <- initial_split(full_dataset, 0.75, strata = subfamily)
train <- training(df_split)

train |> count(subfamily) |> print.AsIs()

rm(full_dataset)


## Setup Nested CV --------------------------------------------------

set.seed(1234)
nest_cv <- 
  nested_cv(
    train, 
    outside = vfold_cv(v = 2, repeats = 2, strata = subfamily),
    inside = vfold_cv(v = 2, repeats = 2, strata = subfamily)
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
    inner_id = glue('{outer_id}_{id}_{id2}'),
    inner_splits = splits
  ) |> 
  nest(inner_resamples = c(inner_id, inner_splits))

nest_cv |> unnest(inner_resamples)

rm(train, df_split)


# MAIN / TRAIN & TEST ----------------------------------------

# **evaluate models with nested CV**
tic()
nest_cv_results <- nest_cv |> fit_nested_cv()
toc()
beepr::beep()

write_rds(nest_cv_results, glue('./results/nest_cv_results.rds'))

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
final_summary <- 
  bind_rows(nest_cv_summary, thresh_res) |> 
  filter(.metric %in% c('mcc')) 

# compare models by MCC
print(final_summary)

library(gt)

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
  tab_header(title = "Nested 2x 2-fold cross-validation results") 




