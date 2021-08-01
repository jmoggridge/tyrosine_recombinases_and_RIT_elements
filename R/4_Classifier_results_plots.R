
library(tidyverse)
library(ggrepel)
library(glue)
library(gt)
library(ggtext)


# 0 DATASET PlOTS ------------------------------------------------------------

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
  # bind_rows(smart, non_int) |> mutate(ds = 'All data'),
  train |> mutate(ds = 'Training'),
  test |> mutate(ds = 'Testing')
) |> 
  filter(subfamily != 'non_integrase') |> 
  select(ds, subfamily) |> 
  count(ds, subfamily) |> 
  ggplot(aes(n, fct_rev(subfamily), fill = ds)) +
  geom_col() +
  theme_bw() +
  theme(legend.position = c(0.8, 0.8)) +
  labs(x = 'n sequences', y = NULL, fill = NULL,
       title = 'Class sizes in filtered dataset')


# combine and label dataset
combo <- bind_rows(
  train |> mutate(ds = 'Training'),
  test |> mutate(ds = 'Testing')
) 

combo |> 
  filter(subfamily == 'Other') |> 
  select(acc, ds) |> 
  left_join(non_int |> select(acc, group)) |> 
  count(ds, group) |> 
  ggplot(aes(x = n, y = fct_rev(group), fill = ds)) +
  geom_col() +
  theme_bw() +
  labs(x = 'n sequences', y = NULL,
       title = 'Subgroup sizes in the non-integrases dataset')


## Show sequence lengths of test and train datasets
combo |> 
  filter(subfamily == 'Other') |> 
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
  filter(subfamily != 'Other') |> 
  mutate(len = nchar(prot_seq)) |> 
  ggplot(aes(len, fct_rev(subfamily))) +
  geom_jitter(size = 0.05, shape = 1, alpha = 0.2) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), size = 0.67, alpha = 0.5,
              trim = T, scale = 'width', fill = NA,
              colour = 'red') +
  theme_bw() +
  labs(x = 'n residues', y = NULL,
       title = 'Sequence lengths in the SMART dataset')



rm(combo, smart, non_int, train, test)

# 1 NESTED CV PLOTS ----------------------------------------------------------

# TODO stuff to take from nested cv
# check out predictions and see which classes are misclassified most often
# sample of HMM scores
## check out model parameters

run_name <- '3x3-fold_07-08'

models <- read_rds('./data/unfitted_parsnip_model_set.rds')
glimpse(models)

nest_cv_results <- read_rds('./results/3x3-fold_07-08_nest_cv_results.rds') |> 
  select(-test, -train)
glimpse(nest_cv_results)


nest_cv_summary <- read_rds('./results/3x3-fold_07-08_nest_cv_summary.rds') |> 
  filter(!model_type == 'decision tree')
glimpse(nest_cv_summary)

#### models selected by inner CV

# best tunings
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

#### table 1 -------

# compare models by metrics
nest_cv_summary |>
  select(-values, -n_folds) |> 
  mutate(model_type = ifelse(model_type == 'multinomial regression',
                             'logistic regression', model_type)) |> 
  mutate(across(.cols = c(mean, err), .fns = ~round(x = .x, digits = 5))) |> 
  transmute(
    model = str_replace_all(model_type, '_', ' '),
    metric = str_replace_all(.metric, '_', ' '),
    mean_sd = glue('{mean} ± {err}')
    ) |> 
  filter(metric %in% c('mcc', 'bal accuracy', 'f meas', 'sens', 'spec')) |> 
  
  pivot_wider(id_cols = model, names_from = metric, values_from = mean_sd) |> 
  gt() |> 
  tab_header(title = "Performance estimate from nested 3-fold cross-validation repeated 3 times") 


#### Inner CV results for plots ----------------------------------------

inner_cv_mcc <- 
  nest_cv_results |> 
  select(outer_id, inner_cv) |> 
  unnest(inner_cv) |>
  filter(.metric == 'mcc') |> 
  select(-.metric) |> 
  left_join(models, by = c("model_type", "model_id"))

inner_cv_mcc



#### knn  --------------------------------------------------------

knn_pooled_scores <- 
  inner_cv_mcc |> 
  filter(str_detect(model_type, 'nearest')) |> 
  unnest(params) |> 
  unnest(values) |> 
  group_by(outer_id, neighbors) |> 
  mutate(fold = as_factor(row_number())) |> 
  ungroup()
knn_pooled_scores

grouped_knn_scores <- knn_pooled_scores |> 
  mutate(values = round(values, 5)) |> 
  count(neighbors, values)

knn_scores_plot <- 
  knn_pooled_scores |> 
  ggplot(aes(y = values, x = neighbors)) +
  geom_point(data = grouped_knn_scores, 
             aes(y = values, x = neighbors, size = n),
             shape = 1) +
  geom_smooth(size = 1, alpha = 0.8, color = 'blue', span = 0.3, se = F,
              method = 'loess', formula = 'y ~ x') +
  theme_bw() +
  theme(legend.text = element_text(size = 7))

knn_scores_plot +
  labs(y = 'MCC', size = 'folds', title = 'K-Nearest Neighbor')

knn_plot <- knn_scores_plot +
  labs(y = 'MCC', size = 'folds', 
       title = 'MCC as a function of k for k-Nearest Neighbor',
       subtitle = 'Points show 81 inner folds; line shows Loess fit over all inner folds.'
       )
knn_plot

rm(knn_pooled_scores, grouped_knn_scores, knn_scores_plot)

#### glmnet -------------------------------

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
  theme_bw() 

glmnet_plot +
  labs(title = 'Elastic net logistic regression',
       fill = 'MCC mean', size = 'MCC sd')

glmnet_plot <- glmnet_plot +
  labs(title = 'MCC as a function of penalty and mixture parameters for elastic net',
       subtitle = 'Colors represent means of 81 inner folds, points represent sd',
       fill = 'MCC mean', size = 'MCC sd')

rm(glmnet_pooled_scores)

#### RF ----
rf_scores <- 
  inner_cv_mcc |> 
  filter(str_detect(model_type, 'forest')) |> 
  unnest(params) |> 
  select(-mean, -err, -n_folds) |> 
  unnest(values) |> 
  group_by(outer_id, mtry) |> 
  mutate(inner_fold = as_factor(paste0(outer_id, '_', row_number()))) |> 
  ungroup() 

rf_scores |> print(n=25)

rf_summary <- rf_scores |> 
  group_by(mtry) |> 
  summarize(
    mean = mean(values, na.rm=T), 
    error = sd(values, na.rm = T),
    .groups = 'drop'
  ) 

grouped_rf_scores <- rf_scores |> 
  mutate(values = round(values, 4)) |> 
  count(mtry, values)

rf_plot <- rf_scores |> 
  ggplot(aes(x = mtry, y = values)) +
  geom_point(data = grouped_rf_scores, shape = 1,
             aes(x = mtry, y = values, size = n)) +
  geom_path(data = rf_summary, 
            aes(x = mtry, y = mean), stat = 'smooth',
            color = 'blue1', span = 1.2,
            se = F, size = 1.5, alpha  = 0.7,
  ) +
  theme_bw() 

rf_plot + 
  labs(x = 'mtry', y = 'MCC', size = 'folds', title = 'Random Forest')

rf_plot <- rf_plot +
  labs(x = 'mtry', y = 'MCC', size = 'folds',
       title = "MCC as a function of random forest hyperparameter *mtry*",
       subtitle = 'Points show 81 inner folds; line shows Loess fit over all inner folds.') +
  theme(plot.title = element_markdown())

rm(grouped_rf_scores, rf_summary, rf_scores)


#### NestCV summary ----------------------

## outer folds results
outer_scores_df <- 
  nest_cv_summary |> 
  filter(model_type !=  "decision tree") |> 
  filter(!.metric %in% c('precision', 'recall', 'kap', 'f_meas', 
                         'accuracy', 'npv', 'ppv')) |> 
  mutate(
    .metric = as.character(.metric),
    .metric = as_factor(case_when(
      .metric == 'sens' ~ 'sensitivity',
      .metric == 'spec' ~ 'specificity',
      TRUE ~ .metric
    ))) |> 
  mutate(
    min = map_dbl(values, min),
    max = map_dbl(values, max),
    .metric = fct_relevel(.metric, 'specificity', 'sensitivity', 'mcc'),
    model_type = ifelse(model_type == 'multinomial regression', 'logistic regression', model_type)
  )


# performance metrics across outer folds with models tuned by the inner cv
all_metrics_plot <- outer_scores_df  |> 
  ggplot(aes(
    y = fct_rev(model_type), 
    x = mean, 
    xmax = mean + err, 
    xmin = mean - err
  )) +
  geom_vline(aes(xintercept = 1), alpha = 0.25) +
  ggbeeswarm::geom_quasirandom(
    data = outer_scores_df |> unnest(values),groupOnX = F,
    aes(x = values), shape = 16, alpha = 0.33
  ) +
  geom_pointrange(color = 'red2', alpha = 0.6, size = 1, fatten = 1.2) +
  facet_wrap(~.metric, nrow = 3, scales = 'free_x') +
  labs(y = NULL, 
       x = NULL, 
       title = '3-fold nested cross-validation repeated 3 times',
       subtitle = 
         'Performance estimates for the hyperparameter tuning process.\nPoint-ranges shows mean +/- sd; points show outer-fold performance') +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 10))

all_metrics_plot

save(knn_plot, rf_plot, glmnet_plot, all_metrics_plot, file = './results/3x3-fold-07-08-nested_cv_plots.rda')



rm(all_metrics_plot, best_tunes, outer_scores_df, knn_plot, rf_plot, glmnet_plot,
   nest_cv_results, nest_cv_summary, inner_cv_mcc)


# 2 CV PLOTS ---------------------------------------------------------------

# summary of CV results
cv_res <- read_rds('results/regular_cv_07-10/cv_res.rds') |> 
  filter(!model_type == 'decision tree')

glimpse(cv_res)

# results for each fold
fold_res <- read_rds('results/regular_cv_07-10/fold_res.rds') |> 
  select(outer_id, res) |> 
  unnest(res) |> 
  filter(!model_type == 'decision tree')

glimpse(fold_res)

best_models <- cv_res |>
  filter(.metric == 'mcc') |>
  group_by(model_type) |>
  filter(mean == max(mean, na.rm = T))

best_models |>
  left_join(models) |>
  split(~model_type) |>
  map(~unnest(.x, params))

## REG CV Plots ------

# plot MCC
mcc <- cv_res |>
  filter(.metric == 'mcc') |>
  # filter(mean > 0.9) |>
  mutate(
    model_id = ifelse(model_id < 10, paste0('0', model_id), paste(model_id)),
    model = glue('{model_type}_{model_id}'))
mcc

# TODO fit thresholds


## knn results
cv_res |>
  filter(.metric == 'mcc', model_type == 'nearest neighbor') |>
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

#### elastic net CV ------
glmnet_mod_labels <- 
  models |> 
  filter(str_detect(model_type, 'regression')) |> 
  select(model_id, params) |> 
  unnest(params)

cv_res |>
  filter(.metric == 'mcc', str_detect(model_type, 'regression')) |>
  filter(mean>0.9) |> 
  left_join(models) |>
  select(-c(.metric, n_folds, spec)) |>
  unnest(params) |>
  ggplot(aes(penalty, mixture)) +
  geom_tile(aes(fill = mean), color = 'black') +
  geom_point(aes(size = err)) +
  scale_fill_viridis_c(begin = 0.2) +
  scale_x_log10() +
  theme_bw() 

rm(glmnet_mod_labels)

#### RF CV -----------

# extract and unnest RF scores from regular CV
rf_scores <- cv_res |>
  filter(.metric == 'mcc', model_type == 'random forest') |>
  left_join(models, by = c("model_type", "model_id")) |>
  unnest(params) |>
  select(-c(.metric, n_folds, spec)) |> 
  unnest(values) |>
  select(-mean, -err) |> 
  group_by(model_id) |>
  mutate(fold = row_number()) |> 
  ungroup()

rf_scores

grouped_rf_scores <- rf_scores |> 
  mutate(values = round(values, 4)) |> 
  count(mtry, values)

rf_scores |> 
  ggplot(aes(mtry, values)) +
  geom_point(data = grouped_rf_scores, shape = 1, 
             aes(mtry, values, size = n)) +
  geom_smooth(se = F) +
  theme_bw() +
  labs(x = 'mtry', y = 'MCC', size = 'folds', 
       title = 'Random Forest tuning 3x3-fold CV')

rm(list= ls())


# 3 FINAL VALIDATION ---------------------------------------------------------

final_res <- read_rds('results/final_model_07-11/final_validation_results.rds') |> 
  unnest(final_metrics)
glimpse(final_res)

fitted_models <- read_rds('results/final_model_07-11/fit_models.rds') |> 
  select(-fitted_wkfl, -reg_cv_res)
  
glimpse(fit_models)


## CONFUSION MATRICES ----------------------------------------------------

## k=7 kNN confusion matrix
knn_conf_mat <- 
  fitted_models |> 
  filter(model_type == 'nearest neighbor') |> 
  unnest(params) |> 
  filter(neighbors == 7) |>
  select(preds) |> 
  unnest(preds) |> 
  # add true labels from test dataset
  bind_cols(read_rds('./data/classif_test_set.rds') |> 
              select(subfamily)) |> 
  count(subfamily, .pred_class) 

knn_conf_mat |> 
  ggplot(aes(y = fct_rev(.pred_class), x = subfamily, size = n, label = n)) +
  geom_point(alpha = 0.4, show.legend = F) +
  ggrepel::geom_text_repel(size = 4) +
  scale_size(trans = 'log10') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = -35, hjust = 0)) +
  labs(x = 'Truth', y = 'Estimate', 
       title = '7-nearest neighbors')
  


## glmnet confusion matrix
glmnet_conf_mat <- 
  fitted_models |> 
  filter(model_type == 'multinomial regression') |> 
  unnest(params) |> 
  filter(model_id == 39) |>
  select(preds) |> 
  unnest(preds) |> 
  # add true labels from test dataset
  bind_cols(read_rds('./data/classif_test_set.rds') |> 
              select(subfamily)) |> 
  count(subfamily, .pred_class) 

glmnet_conf_mat |> 
  ggplot(aes(y = fct_rev(.pred_class), x = subfamily, size = n, label = n)) +
  geom_point(alpha = 0.4, show.legend = F) +
  ggrepel::geom_text_repel(size = 4) +
  # scale_size(trans = 'log2') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = -35, hjust = 0)) +
  labs(x = 'Truth', y = 'Estimate', 
       title = 'Logistic regression')



## rf confusion matrix
rf_conf_mat <- 
  fitted_models |> 
  filter(model_type == 'random forest') |> 
  unnest(params) |> 
  select(preds) |> 
  unnest(preds) |> 
  # add true labels from test dataset
  bind_cols(read_rds('./data/classif_test_set.rds') |> 
              select(subfamily)) |> 
  count(subfamily, .pred_class) 

rf_conf_mat |> 
  ggplot(aes(y = fct_rev(.pred_class), x = subfamily, size = n, label = n)) +
  geom_point(alpha = 0.4, show.legend = F) +
  ggrepel::geom_text_repel(size = 4) +
  # scale_size(trans = 'log2') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = -35, hjust = 0)) +
  labs(x = 'Truth', y = 'Estimate', 
       title = 'Random forest')








# 
# ## Perfomance Summary ------------------------------------------------
# ## outer folds results
# df <- 
#   nest_cv_summary |> 
#   filter(model_type !=  "decision tree") |> 
#   filter(!.metric %in% c('precision', 'recall', 'f_meas', 'kap', 'accuracy')) |> 
#   mutate(
#     .metric = as.character(.metric),
#     .metric = as_factor(case_when(
#       .metric == 'sens' ~ 'sensitivity',
#       .metric == 'spec' ~ 'specificity',
#       TRUE ~ .metric
#     ))) |> 
#   mutate(
#     min = map_dbl(values, min),
#     max = map_dbl(values, max),
#     .metric = fct_relevel(.metric, 'specificity', 'sensitivity', 'mcc')
#   )
# 
# # performance metrics across outer folds with models tuned by the inner cv
# all_metrics_plot <- df  |> 
#   ggplot(aes(
#     y = fct_rev(model_type), 
#     x = mean, 
#     xmax = mean + err, 
#     xmin = mean - err
#   )) +
#   geom_vline(aes(xintercept = 1), alpha = 0.25) +
#   # geom_jitter(data = df |> unnest(values),
#   #             aes(x = values), shape = 1, alpha = 0.7) +
#   ggbeeswarm::geom_quasirandom(
#     data = df |> unnest(values),groupOnX = F,
#     aes(x = values), shape = 16, alpha = 0.33
#   ) +
#   # geom_errorbarh(color = 'blue1', alpha = 0.5) +
#   geom_pointrange(color = 'red2', alpha = 0.6, size = 1, fatten = 1.2) +
#   facet_wrap(~.metric, nrow = 3, scales = 'free_x') +
#   labs(y = NULL, 
#        x = NULL, 
#        title = '3-fold nested cross-validation repeated 3 times',
#        subtitle = 'Performance estimates for the hyperparameter tuning process') +
#   theme_bw() +
#   theme(axis.text.x = element_text(size = 8),
#         axis.text.y = element_text(size = 10))
# 
# all_metrics_plot

