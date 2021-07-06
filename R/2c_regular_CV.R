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
# set of models with grid search
source('./R/00_get_model_specs.R')

# make directories for prep files
out_path <- '3x3_regular_CV_07-02'
make_dirs(out_path)
system(glue('mkdir ./results/{out_path}/'))

# same object forom 2b_nested_cv.R
nest_cv <- read_rds('./results/nest_cv_06-30/3x3-fold_06-30_nest_cv_results.rds')

outer_cv <- nest_cv |> 
  select(outer_id, outer_splits)
outer_cv
rm(nest_cv)


## Do CV -----

# prepare scored datasets for each fold
#  training split > align, build HMM, score sequences, model 
tic()
cv_prep <- outer_cv |> 
  mutate(prep = map2(.x = outer_splits, 
                     .y = outer_id, 
                     .f = ~prep_data(.x, .y, out_path = out_path))
         )
beepr::beep()
toc()
write_rds(cv_prep, glue('./results/{out_path}/cv_prep.rds'), compress = 'gz')

cv_prep <- read_rds(glue('./results/{out_path}/cv_prep.rds'))
# unnest train and test columns
cv_prep <- cv_prep |> 
  unnest(prep) |> 
  select(-outer_splits)

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
write_rds(fold_res, glue('./results/{out_path}/fold_res.rds'), compress = 'gz')


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

write_rds(cv_res, glue('./results/{out_path}/cv_res.rds'), compress = 'gz')

## select which model is best... check plots first too


best_mods <- cv_res |> 
  group_by(model_type) |> 
  filter(.metric == 'mcc') |> 
  filter(mean - err == max(mean - err, na.rm = T)) 
  
write_rds(top_dogs, './results/3x3_regular_CV_07-02/best_models.rds')


## t-test for best models?
t_test_df <- 
  best_mods |> 
  ungroup() |> 
  select(model_type, model_id, values) |> 
  transmute(model = paste0(model_type, model_id),
            values)

t_test_df |>
  set_names(c('model_1', 'values_1')) |> 
  crossing(t_test_df |> set_names(c('model_2', 'values_2')), 
           .name_repair = 'universal') |>  
  filter(model_1 != model_2) |> 
  mutate(comparision = paste0(sort(c(model_1, model_2))))
