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

# same object forom 2b_nested_cv.R
nest_cv <- read_rds('./results/3x3-fold_06-30_nest_cv_results.rds')

outer_cv <- nest_cv |> 
  select(outer_id, outer_splits)

rm(nest_cv)

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

# evaluate model set on each fold
cat('\n', blue$bold("Evaluating models"))
tic()
fold_res <- cv_prep |>
  unnest(prep) |> 
  mutate(rs = map2(train, test, ~ eval_model_set(models, .x, .y))) |>
  unnest(rs)
toc()
beepr::beep()

# summarize metrics across folds
cat('\n', blue$bold("Summarizing results"))

cv_res <- fold_res |>
  group_by(model_type, model_id, .metric) |>
  summarize(
    mean = mean(.estimate, na.rm = T),
    err = sd(.estimate, na.rm = T),
    n_folds = sum(!is.na(.estimate)),
    values = list(.estimate),
    .groups = 'drop'
  )


# plot MCC
cv_res |> 
  filter(.metric == 'mcc') |> 
  ggplot(aes(x = factor(model_id), y = mean, ymin = mean-err, ymax = mean+err)) +
  geom_pointrange() +
  facet_wrap(~model_type, nrow = 3)


best_models <- cv_res |> 
  filter(.metric == 'mcc') |> 
  group_by(model_type) |> 
  filter(mean == max(mean, na.rm = T))

best_models |> 
  left_join(models) |> 
  split(~model_type) |> 
  map(~unnest(.x, params))

# fit thresholds


## knn results
cv_res |> 
  filter(.metric == 'mcc', model_type == 'nearest_neighbor') |> 
  left_join(models) |> 
  unnest(params) |> 
  select(-c(.metric, n_folds, spec)) |> 
  unnest(values) |> 
  group_by(model_id) |> 
  mutate(fold = row_number()) |> 
  ggplot(aes(neighbors, values)) +
  geom_jitter(shape = 1, width = 0.05, height = 0) +
  geom_line(aes(group = fold), stat = 'smooth', se = F, color = 'darkgray', span = 0.3) +
  geom_smooth(se = F) +
  theme_classic()

# elastic net res
cv_res |> 
  filter(.metric == 'mcc', str_detect(model_type, 'glmnet')) |> 
  left_join(models) |> 
  select(-c(.metric, n_folds, spec)) |> 
  unnest(params) |> 
  unnest(values) |> 
  View()
  