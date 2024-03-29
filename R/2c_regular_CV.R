## 2c Regular cross-validation:  3-fold CV, repeated 3 times.


## Setup ----

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

# data prep and modelling functions
source('./R/00_functions.R')

# make directories for prep files
run_name <- 'regular_cv_07-10'
out_path <- glue('{here::here()}/results/{run_name}')
out_path
system(glue('mkdir {out_path}'))
make_dirs(out_path)

# set of models with grid search
models <- read_rds('./data/unfitted_parsnip_model_set.rds')
models

## Data -----

# Outer CV of the same object with splits that was used for nested cv in 2b_nested_cv.R
outer_cv <- 
  read_rds('./results/nested_cv_splits.rds') |> 
  select(outer_id, outer_splits)
outer_cv

## Do CV -----

# prepare scored datasets for each fold
#  training split > align, build HMM, score sequences, model 
tic()
cv_prep <- outer_cv |> 
  mutate(prep = map2(.x = outer_splits, 
                     .y = outer_id, 
                     .f = ~prep_data(.x, .y, out_path = out_path))
         ) |> 
  unnest(prep) |> 
  select(-outer_splits)

# unnest train and test columns
beepr::beep()
toc()

write_rds(cv_prep, glue('{out_path}/cv_prep.rds'), compress = 'gz')

cv_prep <- read_rds(glue('{out_path}/cv_prep.rds'))
cv_prep

# evaluate model set on each fold
cat('\n', blue$bold("Evaluating models"))
tic()
fold_res <- cv_prep |>
  mutate(res = map2(train, test, ~ eval_model_set(models, .x, .y))) 
toc()
beepr::beep()

fold_res
rm(cv_prep)
write_rds(fold_res, glue('{out_path}/fold_res.rds'), compress = 'gz')


# summarize metrics across folds for ML models
cat('\n', blue$bold("Summarizing results"))

cv_res <- fold_res |>
  select(-test, -train) |> 
  unnest(res) |> 
  group_by(model_type, model_id, .metric) |>
  summarize(
    mean = mean(.estimate, na.rm = T),
    err = sd(.estimate, na.rm = T),
    n_folds = sum(!is.na(.estimate)),
    values = list(.estimate),
    .groups = 'drop'
  )
cv_res

# evaluate simple rule-based classification
thresh_res <- fold_res |> 
  mutate(threshold_res = map2(train, test, ~eval_threshold_classififer(.x, .y))) |> 
  select(outer_id, threshold_res) |> 
  unnest(threshold_res) |> 
  mutate(model_type = 'score & threshold rule') |> 
  group_by(model_type, .metric) |>
  summarize(
    mean = mean(.estimate, na.rm = T),
    err = sd(.estimate, na.rm = T),
    n_folds = sum(!is.na(.estimate)),
    values = list(.estimate),
    .groups = 'drop'
  )

cv_res <- bind_rows(cv_res, thresh_res)

write_rds(cv_res, glue('{out_path}/cv_res.rds'), compress = 'gz')

## select which model is best... check plots first too

best_mods <- cv_res |> 
  filter(!str_detect(model_type, 'decision')) |> 
  group_by(model_type) |> 
  filter(.metric == 'mcc') |> 
  filter(mean - err == max(mean - err, na.rm = T))
best_mods
write_rds(best_mods, glue('{out_path}/best_models.rds'))


## TODO t-test to check if differences are significant??
# 
# ## t-test for best models?
# t_test_df <- 
#   best_mods |> 
#   ungroup() |> 
#   select(model_type, model_id, values) |> 
#   transmute(model = glue('{model_type} {model_id}'),
#             values)
# 
#  ttests <- t_test_df |>
#   crossing(t_test_df |> 
#              set_names(c('model_2', 'values_2')),
#            .name_repair = 'minimal') |> 
#   filter(model != model_2)
# 
#  ttests |> unnest(c(values, values_2)) |> print(n=500)
#  
# ttests |> 
#   mutate(
#     test = map2(.x = values, .y = values_2, .f = ~t.test(.x, .y, paired = T)),
#     pval = map_dbl(test, 'p.value')
#     ) |> 
#   arrange(pval)
# 
# 
# t.test(rep(1, 10), rep(1, 10), paired=T)  
# 
