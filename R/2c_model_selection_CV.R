# CV to tune classifiers

# future::plan(multisession, ...) issue
# see issue 511 on future github.
# options(parallelly.makeNodePSOCK.setup_strategy = "sequential")

library(tidyverse)
library(tidymodels)
library(Biostrings)
library(future)
library(furrr)
library(DECIPHER)
library(glue)
library(here)
library(crayon)
library(tictoc)
library(beepr)

source('./code/00_functions.R')

## Directories ----------------------------------------------------------------

# name for directory in project folder to store nested CV files
dir <- 'classifier_02/'

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


## Combine datasets  ------

# integrase dataset
smart_df <- 
  read_rds('./data/SMART/smart_df.rds') |>
  select(subfamily, acc, description, prot_seq, dom_seq) |> 
  mutate(subfamily = as_factor(subfamily))

# negative examples (non-integrases)
non_integrases <- 
  read_rds('./data/non_integrase_seqs/nonint_df.rds') |> 
  mutate(dom_seq = NA,
         subfamily = as_factor(subfamily)) |> 
  select(subfamily, acc, description, prot_seq, dom_seq)

set.seed(123)
# combine datasets 
# TODO don't call this full_Dataset when it's down-sampled! for clarity
combined_dataset <- 
  bind_rows(smart_df, non_integrases) |> 
  # downsample to 10000
  group_by(subfamily) |> 
  slice_sample(n = 15000, replace = F) |> 
  ungroup()

rm(smart_df, non_integrases)

## Initial split -------------------------------------------------

set.seed(54321)

df_split <- initial_split(combined_dataset, 0.75, strata = subfamily)
train <- training(df_split)
train |> count(subfamily) |> print(n=22)

rm(combined_dataset)

## 2x2-fold CV
cv_splits <- 
  vfold_cv(train, v = 2, repeats = 2, strata = subfamily) |> 
  transmute(splits, id = glue('{id}_{id2}'))

cv_splits

fit_cv <- function(cv_splits){
  # preprocess: align, hmms, score seqs
  cat('\n', blue$bold("Preparing alignments, HMMs, & scoring sequences"))
  tic()
  cv_splits <- cv_splits |>
    mutate(prep = map2(splits, id, ~prep_data(.x, .y))) |>
    unnest(prep)

  # fit & evaluate models
  cat('\n', blue$bold("Evaluating models"))
  tic()
  # evaluate model set on each fold
  fold_res <- cv_splits |>
    mutate(rs = map2(train, test, ~eval_model_set(models, .x, .y))) |>
    unnest(rs)
  toc()

  cat('\n', blue$bold("Summarizing results"))
  # summarize metrics across folds
  cv_res <- fold_res |>
    group_by(model_type, model_id, .metric) |>
    summarize(
      mean = mean(.estimate, na.rm = T),
      error = sd(.estimate, na.rm = T),
      n_folds = n(),
      values = list(.estimate),
      .groups = 'drop'
    )
  cv_res
}

cv_res <- fit_cv(cv_splits)


beepr::beep()
