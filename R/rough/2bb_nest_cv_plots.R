library(tidyverse)
library(glue)
library(gt)

source('./R/00_get_model_specs.R')

nest_cv_results <- read_rds('./results/nest_cv_06-30/3x3-fold_06-30_nest_cv_results.rds')
final_summary <- read_rds('./results/nest_cv_06-30/3x3-fold_06-30_nest_cv_summary.rds')

# compare models by metrics
final_summary |>
  select(-values, -n_folds) |> 
  mutate(mean = round(mean, 5),
         err = round(err, 5)) |> 
  transmute(
    model = str_replace_all(model_type, '_', ' '),
    metric = str_replace_all(.metric, '_', ' '),
    mean_sd = glue('{mean} ± \n{err}')) |> 
  filter(metric %in% c('mcc', 'bal_accuracy', 'f_meas', 'kap',
                       'precision', 'recall','sens', 'spec')) |> 
  pivot_wider(id_cols = model, names_from = metric, values_from = mean_sd) |> 
  gt() |> 
  tab_options(table.font.size = 14) |> 
  tab_header(title = "Performance metrics from nested 3-fold cross-validation repeated 3 times") 


## check out model parameters

best_tunes <- nest_cv_results |> 
  select(outer_id, best_tune) |> 
  unnest(best_tune) |> 
  select(-spec) |> 
  left_join(models, by = c("model_type", "model_id")) |> 
  split(~model_type) |> 
  map(~unnest(.x, params)) 

best_tunes

best_tunes |> 
  map(~select(.x, model_type, model_id)) |> 
  bind_rows() |> 
  ggplot(aes(x = model_id)) +
  geom_bar() +
  facet_wrap(~model_type, nrow = 3) +
  theme_bw() +
  scale_y_continuous(breaks = scales::pretty_breaks()) +
  scale_x_continuous(breaks = 1:20, limits = c(0,20)) +
  labs(title = 'Which model tunings were selected most often?',
       x = 'model id', y = 'number of folds') +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor.x = element_blank())

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

# best_scores_on_fold <- knn_pooled_scores |>
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
      title = 'MCC for nearest neighbors classifier in nested CV',
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

glmnet_mod_labels <- 
  models |> 
  filter(str_detect(model_type, 'glmnet')) |> 
  select(model_id, params) |> 
  unnest(params)
  
glment_pooled_scores <- inner_cv_mcc |> 
  filter(str_detect(model_type, 'glmnet')) |> 
  unnest(params) |> 
  unnest(values) |> 
  group_by(outer_id, mixture, penalty) |> 
  mutate(fold = as_factor(row_number())) |> 
  ungroup()

glment_pooled_scores |> 
  group_by(mixture, penalty) |> 
  summarize(mean = mean(values, na.rm=T), .groups = 'drop') |> 
  filter(!is.nan(mean)) |> 
  ggplot(aes(penalty, mixture)) +
  geom_point(aes(size = mean, color = mean)) +
  ggrepel::geom_text_repel(data = glmnet_mod_labels, 
            aes(penalty, mixture, label = model_id)) +
  scale_x_log10() +
  theme_bw() +
  labs(title = 'MCC of Elastic net classifiers in nested CV',
       subtitle = 'Each point represents the mean performance across 81 inner folds',
       color = 'mean MCC', size = NULL)


glment_pooled_scores |> 
  select(model_id, values, fold) |> 
  ggplot(aes(x = as_factor(model_id), y = values)) +
  geom_jitter(alpha = 0.5, shape = 1) +
  geom_boxplot(aes(group = model_id), 
               color = 'red', fill = NA,
               outlier.shape = NA) +  
  labs(x = 'model id', y = 'MCC',
       title = 'Performance of Elastic net classifiers in nested CV',
       subtitle = 'Data for 81 inner folds') +
  theme_bw()


glment_pooled_scores |> 
  select(model_id, values, fold) |> 
  filter(!model_id %in% c(6, 9, 17)) |> 
  ggplot(aes(x = model_id, y = values)) +
  geom_jitter(alpha = 0.5, shape = 1) +
  geom_boxplot(aes(group = model_id), 
               color = 'red', fill = NA,
               outlier.shape = NA) +
  scale_x_continuous(breaks = 1:20, limits = c(0,20)) +
  labs(x = 'model id', y = 'MCC',
       title = 'Performance of Elastic net classifiers in inner CV',
       subtitle = 'Data for 81 inner folds. Showing only models with mean MCC > 0.9.') +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.y = element_blank())





## outer folds results
df <- 
  final_summary |> 
  mutate(
    model = str_replace_all(model_type, '_', ' ') |> str_remove('rpart| glmnet') 
  ) |> 
  filter(!.metric %in% c('npv', 'ppv')) |> 
  mutate(min = map_dbl(values, min),
         max = map_dbl(values, max))

df  |> 
  ggplot(aes(y = model, x = mean, xmax = mean + err, xmin = mean - err)) +
  geom_jitter(data = df |> unnest(values),
              aes(x = values), shape = 1, alpha = 0.75) +
  geom_point(color = 'red') +
  geom_errorbarh(color = 'red') +
  facet_wrap(~.metric, nrow = 3, scales = 'free_x') +
  labs(y = NULL, x = NULL,
       title = 'Performance estimates from applying tuned models in nested CV',
       subtitle = 'Data from 9 outer folds. Grey: individual folds; Red: mean and sd') +

  theme_bw()



# check out predictions and see which classes are misclassified most often


# sample of HMM scores
