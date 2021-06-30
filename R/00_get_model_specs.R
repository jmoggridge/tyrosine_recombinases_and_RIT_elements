library(tidyverse)
library(tidymodels)

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
    grid = grid_max_entropy(
      size = 20, tree_depth(), cost_complexity(), min_n(range = c(2L, 40L))
    ))

# logistic reg models
log_reg_models <- 
  multinom_reg(mode = 'classification') |> 
  set_engine('glmnet') |> 
  make_models(name = 'multinom_reg_glmnet',
              grid = grid_max_entropy(mixture(), penalty(), size = 20))
knn_models <- 
  nearest_neighbor(mode = 'classification') |> 
  set_engine('kknn') |> 
  make_models(name = 'nearest_neighbor',
              grid = grid_regular(neighbors(), levels = 20)
  )
# combine all model specs
models <- bind_rows(tree_models, log_reg_models, knn_models)

rm(tree_models, log_reg_models, knn_models)

