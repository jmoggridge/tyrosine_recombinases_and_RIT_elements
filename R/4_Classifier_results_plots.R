
library(tidyverse)
library(ggrepel)
library(glue)
library(gt)
library(ggtext)
library(tidymodels)
library(themis)
library(patchwork)


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
  geom_boxplot(color = 'red', fill = NA, outlier.alpha = 0) +
  # geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), size = 0.67, alpha = 0.5,
  #             trim = T, scale = 'width', fill = NA,
  #             colour = 'red') +
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

models <- read_rds('./data/unfitted_parsnip_model_set.rds') |> 
  filter(!model_type == 'decision tree')
glimpse(models)

nest_cv_results <- read_rds('./results/3x3-fold_07-08_nest_cv_results.rds') |> 
  select(-test, -train)
glimpse(nest_cv_results)

nest_cv_summary <- read_rds('./results/3x3-fold_07-08_nest_cv_summary.rds') |> 
  filter(!model_type == 'decision tree')
glimpse(nest_cv_summary)

#### models selected by inner CV -----
best_tunes <- nest_cv_results |> 
  select(outer_id, best_tune) |> 
  unnest(best_tune) |> 
  left_join(models, by = c("model_type", "model_id", "spec")) |> 
  split(~model_type) |> 
  map(~unnest(.x, params)) 

best_tunes

#### table 1 -------

# compare models by performance metrics for Outer CV
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


#### NestCV summary ----------------------

## outer folds results
outer_scores_df <- 
  nest_cv_summary |> 
  filter(model_type !=  "decision tree") |> 
  filter(!.metric %in% c('precision', 'recall', 
                         'accuracy', 'npv', 'ppv')) |> 
  mutate(
    .metric = as.character(.metric),
    .metric = as_factor(case_when(
      .metric == 'sens' ~ 'Sensitivity',
      .metric == 'spec' ~ 'Specificity',
      .metric == 'mcc' ~ 'MCC',
      .metric == 'f_meas' ~ 'F Measure',
      .metric == 'kap' ~ "Cohen's Kappa",
      .metric == 'bal_accuracy' ~ 'Balanced Accuracy',
      TRUE ~ .metric
    )),
    model_type = ifelse(model_type == 'multinomial regression', 'logistic regression', model_type),
    model_type = ifelse(model_type == 'maxscore + threshold', 'top-score & threshold', model_type),
    model_type = str_to_title(model_type),
    min = map_dbl(values, min),
    max = map_dbl(values, max),
    .metric = fct_relevel(.metric, 'Specificity', 'Sensitivity', 'MCC'),
  )


# performance metrics across outer folds with models tuned by the inner cv
all_metrics_plot <- outer_scores_df  |> 
  # mutate() |> 
  ggplot(aes(
    y = fct_rev(model_type), 
    x = mean, 
    xmax = mean + err, 
    xmin = mean - err
  )) +
  geom_vline(aes(xintercept = 1), alpha = 0.25) +
  ggbeeswarm::geom_quasirandom(
    data = outer_scores_df |> unnest(values), 
    aes(x = values), shape = 1, alpha = 0.5,
    groupOnX = F
  ) +
  geom_pointrange(color = 'red2', alpha = 0.6, size = 1, fatten = 1.2) +
  facet_wrap(~.metric, nrow = 2, scales = 'free_x') +
  labs(y = NULL, 
       x = NULL, 
       # title = '3-fold nested cross-validation repeated 3 times',
       # subtitle = 
       # 'Performance estimates for the hyperparameter tuning process.\nPoint-ranges shows mean +/- sd; points show outer-fold performance'
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 12), 
        strip.text = element_text(size = 12)
  )

all_metrics_plot

## Statistical testing -----------------------------------------------------

## significance tests for nested CV
MCC_aov <- aov(values ~ model_type, 
      data = 
        outer_scores_df |> 
        unnest(values) |> 
        filter(.metric == 'MCC')
      )
summary(MCC_aov)
# check that variance is homogenous (yes)
# plot(MCC_aov)
TukeyHSD(MCC_aov)



spec_aov <- aov(values ~ model_type, 
               data = 
                 outer_scores_df |> 
                 unnest(values) |> 
                 filter(.metric == 'Specificity')
)
# plot(spec_aov)
summary(spec_aov)
TukeyHSD(spec_aov)

sens_aov <- aov(values ~ model_type, 
                data = 
                  outer_scores_df |> 
                  unnest(values) |> 
                  filter(.metric == 'Sensitivity')
)
# ANOVA test assumes that, the data are normally distributed and the variance across groups are homogeneous:
# variance of sensitivity residuals doesn't conform to assumptions
# ie. there seems to be a relationship between residuals and fitted values
# plot(sens_aov)
summary(sens_aov)
TukeyHSD(sens_aov)

# test for homogeneity of variance
library(car)
leveneTest(values ~ model_type, 
           data = 
             outer_scores_df |> 
             unnest(values) |> 
             filter(.metric == 'Sensitivity')
           )
# p < 0.01; variance across groups is signif different
# use Welch's one-way test instead of anova:
oneway.test(values ~ model_type,
            data = outer_scores_df |>
              unnest(values) |>
              filter(.metric == 'Sensitivity'))
# differences between groups are significant

kruskal.test(values ~ model_type,
            data = outer_scores_df |>
              unnest(values) |>
              filter(.metric == 'Sensitivity'))
# differences between groups are significant


bal_acc_aov <- aov(values ~ model_type, 
                data = 
                  outer_scores_df |> 
                  unnest(values) |> 
                  filter(.metric == 'Balanced Accuracy')
)
summary(bal_acc_aov)
TukeyHSD(bal_acc_aov)

fmeas_aov <- aov(values ~ model_type, 
                   data = 
                     outer_scores_df |> 
                     unnest(values) |> 
                     filter(.metric == 'F Measure')
)
summary(fmeas_aov)
TukeyHSD(fmeas_aov)



rm(sens_aov, spec_aov, MCC_aov, fmeas_aov, bal_acc_aov,
   sens_owt)



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
  theme_bw() + 
  labs(x = 'mtry', y = 'MCC', size = 'folds', title = 'Random Forest')

rf_plot <- rf_plot +
  labs(x = 'mtry', y = 'MCC', size = 'folds',
       title = "MCC as a function of random forest hyperparameter *mtry*",
       subtitle = 'Points show 81 inner folds; line shows Loess fit over all inner folds.') +
  theme(plot.title = element_markdown())


glmnet_plot <- glmnet_plot + labs(title = NULL, subtitle = NULL)
knn_plot <- knn_plot + labs(title = NULL, subtitle = NULL)
rf_plot <- rf_plot + labs(title = NULL, subtitle = NULL)
glmnet_plot + knn_plot + rf_plot + plot_annotation(tag_levels = 'A')


rm(knn_pooled_scores, grouped_knn_scores, knn_scores_plot)
rm(glmnet_pooled_scores)
rm(grouped_rf_scores, rf_summary, rf_scores)




#### end of Nested CV ------

save(knn_plot, rf_plot, glmnet_plot, all_metrics_plot, file = './results/3x3-fold-07-08-nested_cv_plots.rda')

rm(all_metrics_plot, best_tunes, outer_scores_df, knn_plot, rf_plot, glmnet_plot,
   nest_cv_results, nest_cv_summary, inner_cv_mcc)


# 2 CV PLOTS ---------------------------------------------------------------

models <- read_rds('data/unfitted_parsnip_model_set.rds')
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

# best models by MCC
best_models <- cv_res |>
  filter(.metric == 'mcc') |>
  group_by(model_type) |>
  filter(mean == max(mean, na.rm = T))
best_models

# select model params
best_models |>
  full_join(models) |>
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


knn_pooled_scores <- 
  cv_res |>
  filter(.metric == 'mcc', model_type == 'nearest neighbor') |>
  left_join(models) |>
  unnest(params) |>
  select(-c(.metric, n_folds, spec)) |>
  unnest(values) 
knn_pooled_scores

grouped_knn_scores <- knn_pooled_scores |> 
  mutate(values = round(values, 5)) |> 
  count(neighbors, values)

knn_plot <- 
  knn_pooled_scores |> 
  ggplot(aes(y = values, x = neighbors)) +
  geom_point(data = grouped_knn_scores, 
             aes(y = values, x = neighbors, size = n),
             shape = 1) +
  geom_smooth(size = 1, alpha = 0.8, color = 'blue', span = 0.3, se = F,
              method = 'loess', formula = 'y ~ x') +
  theme_bw() +
  labs(y = 'MCC', size = 'folds', title = 'K-Nearest Neighbor')

knn_plot 


#### elastic net CV ------
glmnet_mod_labels <- 
  models |> 
  filter(str_detect(model_type, 'regression')) |> 
  select(model_id, params) |> 
  unnest(params)

glmnet_plot <- cv_res |>
  filter(.metric == 'mcc', str_detect(model_type, 'regression')) |>
  filter(mean>0.999) |> 
  left_join(models) |>
  select(-c(.metric, n_folds, spec)) |>
  unnest(params) |>
  ggplot(aes(penalty, mixture)) +
  geom_tile(aes(fill = mean), color = 'black') +
  geom_point(aes(size = err)) +
  scale_fill_viridis_c(begin = 0.2) +
  scale_x_log10() +
  theme_classic() +
  labs(
    fill = 'Mean MCC',
    size = 'MCC sd',
    title = 'Logistic Regression')

glmnet_plot

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

rf_plot <- rf_scores |> 
  ggplot(aes(mtry, values)) +
  geom_point(data = grouped_rf_scores, shape = 1, 
             aes(mtry, values, size = n)) +
  geom_smooth(se = F) +
  theme_bw() +
  labs(x = 'mtry', y = 'MCC', size = 'folds', 
       title = 'Random Forest')



glmnet_plot <- glmnet_plot + labs(title = NULL) 
knn_plot <- knn_plot + labs(title = NULL) 
rf_plot <- rf_plot + labs(title = NULL) 

glmnet_plot + knn_plot + rf_plot

rm(list= ls())


# 3 FINAL VALIDATION ---------------------------------------------------------

final_res <- read_rds('results/final_model_07-11/final_validation_results.rds') |> 
  unnest(final_metrics)
glimpse(final_res)

fitted_models <- read_rds('results/final_model_07-11/fit_models.rds') |> 
  select(-fitted_wkfl, -reg_cv_res)
  
prepped_data <- read_rds('./results/final_model_07-11/prepped_data.rds') 

## Threshold rule classifier -----
train_prep <- prepped_data |>  
  pull(train_prep) |>
  pluck(1) |> 
  select(-contains('seq'))
test_prep <-prepped_data |>  
  pull(test_prep) |>
  pluck(1) |> 
  select(-contains('seq'))

# data preprocessing
# specify modelling workflow with scaling and ignore seqs and id
recip <- 
  recipe(subfamily ~ ., data = train_prep) |> 
  update_role(acc, new_role = "id") |> 
  step_smote(subfamily, over_ratio = 0.25) |> 
  step_normalize(all_predictors())

# normalized + smoted thresholds
normalized_and_smoted <- recip |> prep() |> juice()

# eg. a training set with SMOTE done
normalized_and_smoted |> 
  count(subfamily) |> 
  ggplot(aes(n, subfamily)) + geom_col() +
  labs(title = 'Training set with SMOTE') +
  theme_bw()

# set thresholds
normalized_smoted_thresholds <-
  normalized_and_smoted |> 
  select(subfamily, Arch1:Xer) |> 
  pivot_longer(Arch1:Xer, names_to = 'HMM', values_to = 'score') |> 
  filter(HMM == subfamily) |> 
  group_by(HMM) |> 
  filter(score == min(score)) |> 
  transmute(HMM, threshold = score)

normalized_smoted_thresholds
normalized_smoted_thresholds |> ggplot(aes(threshold, HMM)) + geom_col()
write_rds(normalized_smoted_thresholds,
          './models/hmm_normalized_smoted_thresholds.rds')

normalized_and_smoted |> 
  relocate(subfamily) |> 
  filter(round(Int_Tn916,3) > -0.210) |> 
  count(subfamily) |> 
  ggplot(aes(n, subfamily)) +
  geom_col()

source('R/00_functions.R')
thresholds <- set_thresholds(train_prep)
classed_test <- hmm_threshold_class(test = test_prep, thresholds = thresholds)
# collect metrics
class_metrics <- metric_set(yardstick::mcc, kap, sens, spec, 
                            f_meas, ppv, npv, accuracy, bal_accuracy)
classed_test |>
  class_metrics(subfamily, estimate = .pred_class)



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
              select(subfamily, acc)) |> 
  group_by(subfamily, .pred_class) |> 
  summarize(n = length(acc),
            acc_list = list(acc)) |> 
  ungroup() |> 
  mutate(model = 'knn')
  
## glmnet confusion matrix
glmnet_conf_mat <- 
  fitted_models |> 
  filter(model_type == 'multinomial regression') |> 
  unnest(params) |> 
  filter(model_id == 30) |>
  select(preds) |> 
  unnest(preds) |> 
  # add true labels from test dataset
  bind_cols(read_rds('./data/classif_test_set.rds') |> 
              select(subfamily, acc)) |> 
  group_by(subfamily, .pred_class) |> 
  summarize(n = length(acc),
            acc_list = list(acc)) |> 
  ungroup() |> 
  mutate(model = 'logistic regression')

## rf confusion matrix
rf_conf_mat <- 
  fitted_models |> 
  filter(model_type == 'random forest') |> 
  unnest(params) |> 
  select(preds) |> 
  unnest(preds) |> 
  # add true labels from test dataset
  bind_cols(read_rds('./data/classif_test_set.rds') |> 
              select(subfamily, acc)) |> 
  group_by(subfamily, .pred_class) |> 
  summarize(n = length(acc),
            acc_list = list(acc)) |> 
  ungroup() |> 
  mutate(model = 'random forest')

# are the errors all the same sequences?
bind_rows(
  knn = knn_conf_mat |> 
    filter(subfamily != .pred_class) |> 
    unnest(acc_list),
  glmnet = glmnet_conf_mat |> 
    filter(subfamily != .pred_class) |> 
    unnest(acc_list),
  rf = rf_conf_mat |> 
    filter(subfamily != .pred_class) |> 
    unnest(acc_list)) |> 
  arrange(subfamily, .pred_class)

# the xer/integron one is the same
# the in tn916/p2 preds are the same in all
# the 916/BPP is found in knn and glmnet

bind_rows(
  classed_test |> 
    count(subfamily, .pred_class) |> 
    mutate(model = 'Top score & threshold rule'),
  knn_conf_mat |> 
    mutate(model = '7-nearest neighbors'),
  glmnet_conf_mat |> 
    mutate(model = 'Logistic regression (alpha = 0.25, lambda = 1e-5)'),
  rf_conf_mat |> 
    mutate(model = 'Random forest (mtry = 6)'),
  ) |> 
  ggplot(aes(subfamily, fct_rev(.pred_class), size = n, label = n)) +
  geom_point(alpha = 0.2, show.legend = F) +
  ggrepel::geom_text_repel(size = 3) +
  facet_wrap(~model) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = -90, hjust = 0)) +
  labs(x = 'Truth', y = 'Estimate')  

## four confusion matrix plots
conf_rules <- classed_test |> 
  count(subfamily, .pred_class) |> 
  ggplot(aes(subfamily, fct_rev(.pred_class), size = n, label = n)) +
  geom_point(alpha = 0.2, show.legend = F) +
  ggrepel::geom_text_repel(size = 3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = -90, hjust = 0)) +
  labs(x = NULL, y = 'Estimate', 
       subtitle = 'Top score & threshold rules')

conf_knn <- knn_conf_mat |> 
  ggplot(aes(y = fct_rev(.pred_class), x = subfamily, size = n, label = n)) +
  geom_point(alpha = 0.2, show.legend = F) +
  ggrepel::geom_text_repel(size = 3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = -90, hjust = 0)) +
  labs(x = NULL, y = NULL, 
       subtitle = '7-nearest neighbors')

conf_glmnet <- glmnet_conf_mat |> 
  ggplot(aes(y = fct_rev(.pred_class), x = subfamily, size = n, label = n)) +
  geom_point(alpha = 0.2, show.legend = F) +
  ggrepel::geom_text_repel(size = 3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = -90, hjust = 0)) +
  labs(x = 'Truth', y = 'Estimate', 
       subtitle = 'Logistic regression (alpha = 0.25, lambda = 1e-5)')

conf_rf <- rf_conf_mat |> 
  ggplot(aes(y = fct_rev(.pred_class), x = subfamily, size = n, label = n)) +
  geom_point(alpha = 0.2, show.legend = F) +
  ggrepel::geom_text_repel(size = 3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = -90, hjust = 0)) +
  labs(x = 'Truth', y = NULL, 
       subtitle = 'Random forest (mtry = 6)')

conf_rules + conf_knn + conf_glmnet + conf_rf


## final test performance metrics -----------------------

# which models were used?
final_res |> 
  select(model_type, model_id) |> 
  distinct() |> 
  left_join(models) |> 
  unnest(params)

# clean up final test performance metrics
final_rs <- final_res |> 
    filter(model_type !=  "decision tree") |> 
    filter(!.metric %in% c('precision', 'recall', 
                           'accuracy', 'npv', 'ppv')) |> 
  mutate(
    .metric = as.character(.metric),
    .metric = as_factor(case_when(
      .metric == 'sens' ~ 'Sensitivity',
      .metric == 'spec' ~ 'Specificity',
      .metric == 'mcc' ~ 'MCC',
      .metric == 'f_meas' ~ 'F Measure',
      .metric == 'kap' ~ "Cohen's Kappa",
      .metric == 'bal_accuracy' ~ 'Balanced Accuracy',
      TRUE ~ .metric
    )),
    model_type = ifelse(model_type == 'multinomial regression', 'logistic regression', model_type),
    model_type = ifelse(model_type == 'score & threshold rule', 'top-score & threshold', model_type),
    model_type = str_to_title(model_type),
    .metric = fct_relevel(.metric, 'Specificity', 'Sensitivity', 'MCC', 
                          'Balanced Accuracy', 'F Measure', "Cohen's Kappa")
    )

# plot performance metrics for final test vs. nested CV
all_metrics_plot2 <-
  outer_scores_df  |> 
  # mutate() |> 
  ggplot(aes(
    y = fct_rev(model_type), 
    x = mean, 

  )) +
  geom_vline(aes(xintercept = 1), alpha = 0.25) +
  ggbeeswarm::geom_quasirandom(
    data = outer_scores_df |> unnest(values), 
    aes(x = values), shape = 1, alpha = 0.5,
    groupOnX = F
  ) +
  geom_pointrange(
    aes(xmax = mean + err, xmin = mean - err),
    color = 'red2', alpha = 0.6, size = 1, fatten = 1.2
    ) +
  geom_point(
    data = final_rs |> transmute(model_type, .metric, .estimate),
    aes(x = .estimate), color = 'black'
    ) +
  facet_wrap(~.metric, nrow = 2, scales = 'free_x') +
  labs(y = NULL, 
       x = NULL, 
       # title = 'Classifier performance in final testing vs. nested CV',
       # subtitle = 
       # 'Black points show final test performance. \nPoint-ranges shows nested CV mean +/- sd; points show values for each outer fold of nested CV; '
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 12), 
        strip.text = element_text(size = 12)
  )

all_metrics_plot2


# 4 Finished Models ---------------------------------------------------------
# check out glmnet classifier coefficients
# for each class, may use scores from multiple HMMs to make a decision
glmnet_mod <- read_rds('models/glmnet_classifier.rds')
glmnet_mod |> 
  pull_workflow_fit() %>%
  tidy() |> 
  View()




