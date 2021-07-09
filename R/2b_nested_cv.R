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

# combine model summary and rule-based classifier summary
final_summary <- bind_rows(nest_cv_summary, thresh_res) 

write_rds(final_summary, glue('./results/{run_name}_nest_cv_summary.rds'))

final_summary |>
  filter(.metric %in% c('mcc', 'bal_accuracy', 'f_meas', 'kap',
                       'precision', 'recall','sens', 'spec')) |>
  View()
  

library(gt)


# compare models by metrics
final_summary |>
  select(-values, -n_folds) |> 
  mutate(across(.cols = c(mean, err), .fns = ~round(x = .x, digits = 5))) |> 
  filter()
  transmute(
    model = str_replace_all(model_type, '_', ' '),
    metric = str_replace_all(.metric, '_', ' '),
    mean_sd = glue('{mean} ± {err}')) |> 
  filter(metric %in% c('mcc', 'bal_accuracy', 'f_meas', 'kap',
                       'precision', 'recall','sens', 'spec')) |> 
  pivot_wider(id_cols = model, names_from = metric, values_from = mean_sd) |> 
  gt() |> 
  tab_header(title = "Performance metrics from nested 3-fold cross-validation repeated 3 times") 


# TODO stuff to take from nested cv
# check out predictions and see which classes are misclassified most often
# sample of HMM scores
## check out model parameters

best_tunes <- nest_cv_results |> 
  select(outer_id, best_tune) |> 
  unnest(best_tune) |> 
  left_join(models, by = c("model_type", "model_id", "spec")) |> 
  split(~model_type) |> 
  map(~unnest(.x, params)) 


best_tunes |> 
  map(~select(.x, model_type, model_id)) |> 
  bind_rows() |> 
  ggplot(aes(x = model_id)) +
  geom_bar() +
  facet_wrap(~model_type, nrow = 3) +
  theme_bw() +
  scale_y_continuous(breaks = scales::pretty_breaks()) +
  labs(title = 'Which model tunings were selected most often?',
       x = 'model id')
## Plot performance as a function of hyperparameters -----

inner_cv_mcc <- 
  nest_cv_results |> 
  select(outer_id, inner_cv) |> 
  unnest(inner_cv) |>
  filter(.metric == 'mcc') |> 
  select(-.metric) |> 
  left_join(models, by = c("model_type", "model_id"))
inner_cv_mcc |> count(model_type, model_id)
## knn scores
knn_pooled_scores <- 
  inner_cv_mcc |> 
  filter(str_detect(model_type, 'nearest')) |> 
  unnest(params) |> 
  unnest(values) |> 
  group_by(outer_id, neighbors) |> 
  mutate(fold = as_factor(row_number())) |> 
  ungroup()
   
# best_scores_on_fold <- knn_scores |> 
#   group_by(fold) |> 
#   filter(values == max(values)) 
  
knn_folds_scores <- 
  knn_scores |> 
  mutate(fold = glue('{outer_id} fold {fold}')) |> 
  split(~fold) |> 
  map(~select(.x, values, neighbors))
  
library(ggtext)

plot_mcc_by_k <- function(pooled, folds) {
  
  p <- knn_pooled_scores |> 
    ggplot(aes(neighbors, values)) +
    labs(
      y = 'MCC', 
      title = 'MCC as a function of *k* parameter for nearest neighbors classifier',
      subtitle = 'Each inner fold is shown as a grey trace, the pooled data is shown in blue'
      ) +
    theme_bw() +
    theme(plot.title = element_markdown(lineheight = 1.2))
  
  for (i in seq_along(knn_folds_scores))
    p <- p + 
      geom_line(data = knn_folds_scores[[i]], stat = 'smooth', alpha = 0.35, se = F)
  
  p <- p +
    geom_line(stat = 'smooth', size = 1.85, alpha = 0.8, color = 'blue', span = 0.3,
              method = 'loess', formula = 'y ~ x')
  p
}

plot_mcc_by_k(knn_pooled_scores, knn_folds_scores)



## glmnet scores

glmnet_pooled_scores <- inner_cv_mcc |> 
  filter(str_detect(model_type, 'glmnet')) |> 
  unnest(params) |> 
  unnest(values) |> 
  group_by(outer_id, mixture, penalty) |> 
  mutate(fold = as_factor(row_number())) |> 
  ungroup()

glmnet_pooled_scores |> 
  group_by(mixture, penalty) |> 
  summarize(mean = mean(values, na.rm=T), .groups = 'drop') |> 
  filter(!is.nan(mean)) |> 
  ggplot(aes(penalty, mixture)) +
  geom_point(aes(size = mean, color = mean)) +
  scale_x_log10() +
  theme_bw() +
  labs(title = 'Elastic net classifier, maximum entropy grid search',
       subtitle = '',
       color = 'mean MCC', size = 'mean MCC')


glmnet_pooled_scores |> 
  select(model_id, values, fold) |> 
  ggplot(aes(x = as_factor(model_id), y = values)) +
  geom_boxplot()
  # geom_jitter(aes(color = fold))


## outer folds results
df <- 
  final_summary |> 
  mutate(
    model = case_when(
      model_type == 'decision_tree_rpart' ~ 'Decision tree',
      model_type == 'multinom_reg_glmnet' ~ 'Elastic net',
      model_type == 'nearest_neighbor' ~ 'k-Nearest neighbor',
      TRUE ~ 'Score + threshold rules',
    )) |> 
  filter(!.metric %in% c('npv', 'ppv')) |> 
  mutate(min = map_dbl(values, min),
         max = map_dbl(values, max))

# performance metrics across outer folds with models tuned by the inner cv
df  |> 
  ggplot(aes(
    y = fct_rev(model), 
    x = mean, 
    xmax = mean + err, 
    xmin = mean - err
    )) +
  geom_jitter(data = df |> unnest(values),
              aes(x = values), shape = 1, alpha = 0.7) +
  # geom_pointrange(fatten = 0.5) +
  geom_errorbarh() +
  facet_wrap(~.metric, nrow = 3, scales = 'free_x') +
  labs(y = NULL, x = NULL) +
  theme_gray()









