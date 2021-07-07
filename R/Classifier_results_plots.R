library(tidyverse)
library(ggrepel)
library(glue)

smart <- read_rds('./data/SMART/smart_df.rds')
non_int <- read_rds('./data/non_integrase_seqs/nonint_df.rds')
train <- read_rds('./data/classif_train_set.rds')
test <- read_rds('./data/classif_test_set.rds')


## Subfamilies & Sequences -----
# show composition of SMART
smart |> 
  count(subfamily) |> 
  ggplot(aes(n, fct_rev(subfamily))) +
  geom_col() +
  theme_bw() +
  labs(x = 'n sequences', y = NULL,
       title = 'Subfamily sizes in the SMART integrases dataset')

# show composition of non_integrase
non_int |> 
  count(group) |> 
  ggplot(aes(n, fct_rev(group))) +
  geom_col() +
  theme_bw() +
  labs(x = 'n sequences', y = NULL,
       title = 'Subfamily sizes in the non-integrases dataset')

## Show composition of test and train datasets

bind_rows(
  smart |> mutate(ds = 'All data'),
  train |> mutate(ds = 'Training'),
  test |> mutate(ds = 'Testing')
) |> 
  filter(subfamily != 'non_integrase') |> 
  select(ds, subfamily) |> 
  count(ds, subfamily) |> 
  ggplot(aes(n, fct_rev(subfamily))) +
  geom_col() +
  facet_wrap(~ds, scales = 'free_x') +
  theme_bw() +
  labs(x = 'n sequences', y = NULL,
       title = 'Subfamily sizes between assessment and validation datasets')


# combine and label dataset
combo <- bind_rows(
  train |> mutate(ds = 'Training'),
  test |> mutate(ds = 'Testing')
) 

combo |> 
  filter(subfamily == 'non_integrase') |> 
  select(acc, ds) |> 
  left_join(non_int |> select(acc, group)) |> 
  count(ds, group) |> 
  ggplot(aes(n, fct_rev(group))) +
  geom_col() +
  facet_wrap(~ds, scales = 'free_x') +
  theme_bw() +
  labs(x = 'n sequences', y = NULL,
       title = 'Subgroup sizes in the non-integrases dataset')


## Show sequence lengths of test and train datasets
combo |> 
  filter(subfamily == 'non_integrase') |> 
  left_join(non_int |> select(acc, group)) |> 
  mutate(len = nchar(prot_seq)) |> 
  ggplot(aes(len, fct_rev(group))) +
  geom_jitter(size = 0.05, shape = 1, alpha = 0.2) +
  geom_boxplot(colour = 'deeppink', outlier.shape = NA, fill = NA) +
  scale_x_log10() +
  theme_bw() +
  labs(x = 'n residues', y = NULL,
       title = 'Sequence lengths in the non-integrases dataset')
  
combo |> 
  filter(subfamily != 'non_integrase') |> 
  mutate(len = nchar(prot_seq)) |> 
  ggplot(aes(len, fct_rev(subfamily))) +
  geom_jitter(aes(color = ds), size = 0.05, shape = 1, alpha = 0.2) +
  # geom_boxplot(colour = 'deeppink', outlier.shape = NA, fill = NA) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), size = 0.67, alpha = 0.8,
              trim = T, scale = 'width', fill = NA,
              colour = 'black') +
  scale_x_log10() +
  theme_bw() +
  guides(legend = override.aes(list = list(size = 5))) +
  labs(x = 'n residues', y = NULL,
       title = 'Sequence lengths in the SMART dataset')



library(tidyverse)
library(glue)

source('./R/00_get_model_specs.R')
cv_res <- read_rds('./results/3x3_regular_CV_07-02/cv_res.rds')


## NEST CV plots ------
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



## REG CV Plots ------

# plot MCC
mcc <- cv_res |> 
  filter(.metric == 'mcc') |> 
  # filter(mean > 0.9) |> 
  mutate(
    model_id = ifelse(model_id < 10, paste0('0', model_id), paste(model_id)),
    model = glue('{model_type}_{model_id}'))

ggplot(mcc) +
  # geom_jitter(data = mcc |> unnest(values), 
  #            aes(factor(model_id), values), 
  #            shape = 1, alpha = 0.5, size = 0.6) +
  geom_pointrange(aes(y = fct_rev(model), x = mean, xmin = mean-err, xmax = mean+err),
                  fatten = 0.3) +
  theme_bw()

sens <- cv_res |> 
  filter(.metric == 'sens') |> 
  filter(mean > 0.8) |> 
  mutate(
    model_id = ifelse(model_id < 10, paste0('0', model_id), paste(model_id)),
    model = glue('{model_type}_{model_id}'))

ggplot(sens) +
  # geom_jitter(data = mcc |> unnest(values), 
  #            aes(factor(model_id), values), 
  #            shape = 1, alpha = 0.5, size = 0.6) +
  geom_pointrange(aes(y = fct_rev(model), x = mean, xmin = mean-err, xmax = mean+err),
                  fatten = 0.3) +
  theme_bw()



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
  theme_classic() +
  labs(y = 'MCC', title = 'k-Nearest neighbor classifier performance in 3-fold CV repeated 3 times',
       subtitle = 'Gray traces: individual folds; blue: pooled data')

# elastic net res
glmnet_mod_labels <- 
  models |> 
  filter(str_detect(model_type, 'glmnet')) |> 
  select(model_id, params) |> 
  unnest(params)

cv_res |>
  filter(.metric == 'mcc', str_detect(model_type, 'glmnet')) |>
  left_join(models) |>
  select(-c(.metric, n_folds, spec)) |>
  unnest(params) |>
  ggplot(aes(penalty, mixture)) +
  geom_point(aes(size = mean, color = mean)) +
  ggrepel::geom_text_repel(data = glmnet_mod_labels, 
                           aes(penalty, mixture, label = model_id)) +
  scale_x_log10() +
  theme_bw() +
  labs(title = 'MCC of Elastic net classifiers in in 3-fold CV repeated 3 times',
       subtitle = 'Each point represents the mean performance across 81 inner folds',
       color = 'mean MCC', size = NULL)


rm(list= ls())

## Final validation results

final_res <- read_rds('./results/final_validation_results_07-05.rds')
fitted_models <- read_rds(glue('./results/final_model_07-05/fitted_models.rds'))



fitted_models |> 
  select(model_type, model_id, final_metrics) |> 
  unnest(final_metrics) |> 
  filter(.metric == 'mcc') |> 
  ggplot(aes(x = .estimate, y = fct_rev(factor(model_id)), fill = model_type)) +
  facet_wrap(~model_type) +
  geom_col(show.legend = F) +
  labs(y = 'model_id', x = 'MCC')

## k=10 NN confusion matrix
my_model <- 
  fitted_models |> 
  filter(model_type == 'nearest_neighbor') |> 
  unnest(params) |> 
  filter(neighbors == 10) 

confmat <- my_model |>
  select(preds) |> 
  unnest(preds) |> 
  # add true labels from test dataset
  bind_cols(read_rds('./data/classif_test_set.rds') |> 
              select(subfamily)) |> 
  count(subfamily, .pred_class) 

confmat |> 
  ggplot(aes(y = fct_rev(.pred_class), x = subfamily, size = n, label = n)) +
  geom_point(alpha = 0.4, show.legend = F) +
  ggrepel::geom_text_repel(size = 4) +
  scale_size(trans = 'log10') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = -35, hjust = 0)) +
  labs(x = 'Truth', y = 'Estimate', 
       title = '10-nearest neighbors classification of final holdout dataset')
  

# glmnet_final_model_set_params ## AGAIN NOT CONSISTENT!!

# fitted_models |> 
#   filter(str_detect(model_type, 'glmnet')) |> 
#   select(contains('model'), params) |> 
#   unnest(params) |> 
#   ggplot(aes(penalty, mixture)) +
#   geom_point() +
#   geom_text_repel(aes(penalty, mixture, label = model_id)) +
#   scale_x_log10() +
#   theme_bw()

