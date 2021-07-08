library(tidyverse)
library(tidymodels)

# TODO need to make it so that the models don't change if rerun!
# set.seed(123)

options(tibble.print_max = 50, tibble.print_min = 20)

### Models  -------------------------------------------------------------------

# update model specifications using parameters grid
make_models <- function(model, grid, name){
  models <- map(
    seq_len(nrow(grid)), ~{
      mods <- model |> 
        update(grid[.x, ])
    })
  tibble(spec = models) |> 
    mutate(model_type = name, 
           model_id = row_number()) |> 
    bind_cols(grid) |> 
    nest(params = names(grid))
}

# decision tree models
tree_models <- 
  decision_tree(mode = 'classification') |> 
  set_engine('rpart') |> 
  make_models(
    name = 'decision_tree_rpart',
    grid = grid_regular(
      levels = 5, 
      tree_depth(range = c(10,30)), 
      cost_complexity(), 
      min_n(range = c(2L, 40L))
    ))
tree_models |> unnest(params) |> 
  ggplot(aes(tree_depth, cost_complexity))  +
  geom_point() +
  scale_y_log10() +
  facet_grid(~min_n)

# logistic reg models
glmnet_models <- 
  multinom_reg(mode = 'classification') |> 
  set_engine('glmnet') |> 
  make_models(
    name = 'multinom_reg_glmnet',
    grid = grid_regular(
      mixture(), 
      penalty(), 
      levels = 12
    ))
glmnet_models |> unnest(params) |> 
  ggplot(aes(penalty, mixture)) + 
  geom_tile(fill = NA, color = 'black') +
  scale_x_log10()

knn_models <- 
  nearest_neighbor(mode = 'classification') |> 
  set_engine('kknn') |> 
  make_models(
    name = 'nearest_neighbor',
    grid = grid_regular(
      neighbors(range = c(2L, 30L)), 
      levels = 29
    ))
knn_models |> unnest(params) |> ggplot(aes(neighbors, 1)) + geom_point()

# combine all model specs
models <- bind_rows(tree_models, glmnet_models, knn_models)
write_rds(models, './data/unfitted_parsnip_model_set.rds')

rm(list = ls())
