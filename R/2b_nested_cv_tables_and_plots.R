
library(tidyverse)
library(gt)
library(glue)

run_name <- '3x3-fold_07-08'

models <- read_rds('./data/unfitted_parsnip_model_set.rds')
nest_cv_results <- read_rds('./results/3x3-fold_07-08_nest_cv_results.rds')
final_summary <- read_rds('./results/3x3-fold_07-08_nest_cv_summary.rds')

# compare models by metrics
final_summary |>
  select(-values, -n_folds) |> 
  mutate(across(.cols = c(mean, err), .fns = ~round(x = .x, digits = 5))) |> 
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
  ggplot(aes(x = as_factor(model_id))) +
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

inner_cv_mcc

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
  knn_pooled_scores |> 
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


knn_plot <- plot_mcc_by_k(knn_pooled_scores, knn_folds_scores)
knn_plot


## glmnet scores

glmnet_pooled_scores <- 
  inner_cv_mcc |> 
  filter(str_detect(model_type, 'reg')) |> 
  unnest(params) |> 
  select(-mean, -err, -n_folds) |> 
  unnest(values) |> 
  group_by(outer_id, mixture, penalty) |> 
  mutate(fold = as_factor(row_number())) |> 
  ungroup()

glmnet_plot <- 
  glmnet_pooled_scores |> 
  group_by(mixture, penalty) |> 
  summarize(
    mean = mean(values, na.rm=T), 
    error = sd(values, na.rm = T),
    .groups = 'drop'
  ) |> 
  filter(!is.nan(mean)) |> 
  filter(mean > 0.9) |> 
  ggplot(aes(penalty, mixture)) +
  geom_tile(aes(fill = mean), color = 'black') +
  geom_point(aes(size = error)) +
  scale_fill_viridis_c(begin = 0.2) +
  scale_x_log10() +
  theme_bw() +
  labs(title = 'MCC as a function of penalty and mixture parameters for elastic net',
       subtitle = 'Colors represent means of 81 inner folds, points represent sd',
       fill = 'MCC mean', size = 'MCC sd')

glmnet_plot

# final_summary |> pull(model_type) |> unique()

## outer folds results
df <- 
  final_summary |> 
  filter(model_type !=  "decision tree") |> 
  filter(!.metric %in% c('npv', 'ppv')) |> 
  mutate(min = map_dbl(values, min),
         max = map_dbl(values, max))

# performance metrics across outer folds with models tuned by the inner cv
all_metrics_plot <- df  |> 
  ggplot(aes(
    y = fct_rev(model_type), 
    x = mean, 
    xmax = mean + err, 
    xmin = mean - err
  )) +
  geom_vline(aes(xintercept = 1), alpha = 0.25) +
  # geom_jitter(data = df |> unnest(values),
  #             aes(x = values), shape = 1, alpha = 0.7) +
  ggbeeswarm::geom_quasirandom(
    data = df |> unnest(values),groupOnX = F,
    aes(x = values), shape = 16, alpha = 0.33
  ) +
  # geom_errorbarh(color = 'blue1', alpha = 0.5) +
  geom_pointrange(color = 'red2', alpha = 0.6, size = 1, fatten = 1.2) +
  facet_wrap(~.metric, nrow = 3) +
  labs(y = NULL, x = NULL, 
       title = '3-fold nested cross-validation repeated 3 times',
       subtitle = 'Performance estimates for the hyperparameter tuning process') +
  theme_bw()

all_metrics_plot


rf_scores <- 
  inner_cv_mcc |> 
  filter(str_detect(model_type, 'forest')) |> 
  unnest(params) |> 
  select(-mean, -err, -n_folds) |> 
  unnest(values) |> 
  group_by(outer_id, mtry) |> 
  mutate(inner_fold = as_factor(paste0(outer_id, '_', row_number()))) |> 
  ungroup() 
rf_scores |> print(n=100)

rf_summary <- rf_scores |> 
  group_by(mtry) |> 
  summarize(
    mean = mean(values, na.rm=T), 
    error = sd(values, na.rm = T),
    .groups = 'drop'
  ) 


rf_plot <- rf_scores |> 
  ggplot(aes(x = mtry, y = values)) +
  geom_path(aes(group = inner_fold), stat = 'smooth',
            size = 0.5, alpha = 0.25, color = 'black') +
  geom_jitter(height = 15e-6,  width = 0,
              size = 0.7, alpha = 0.5, shape = 1) +
  geom_path(data = rf_summary, 
            aes(x = mtry, y = mean), stat = 'smooth',
            color = 'blue1',
            se = F, size = 1.5, alpha  = 0.7,
            ) +
  theme_bw() +
  labs(x = 'mtry', y = 'MCC', 
       title = "MCC as a function of random forest hyperparameter *mtry*",
       subtitle = 'Gray traces and points represent the results for 81 individual inner CV folds;\nthe blue trace shows the average across all inner folds.') +
  theme(plot.title = element_markdown())
  
rf_plot

save(knn_plot, rf_plot, glmnet_plot, all_metrics_plot, file = './results/3x3-fold-07-08-nested_cv_plots.rda')



