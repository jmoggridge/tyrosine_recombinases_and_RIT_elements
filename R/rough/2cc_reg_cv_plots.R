library(tidyverse)
library(glue)

source('./R/00_get_model_specs.R')
cv_res <- read_rds('./results/3x3_regular_CV_07-02/cv_res.rds')


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



