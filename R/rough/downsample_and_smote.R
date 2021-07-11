# downsample and smote

library(tidyverse)
library(tidymodels)
library(themis)
source('./R/00_functions.R')

read_rds('./unit_test/inner_split_train.rds')
train <- read_rds('./data/classif_train_set.rds') |> 
  group_by(subfamily) |> 
  slice_sample(n=100) |> 
  ungroup()

df <- tibble(y = sample(c('a','b'), 100, replace = T), 
            x1 = rnorm(100), x2 = rexp(100)) |> 
  janitor::clean_names()

# specify formula and scaling recipe
recip1 <- 
  recipe(y ~ ., data = df) |> 
  step_nearmiss(y, under_ratio = 1) |> 
  step_normalize(all_predictors())

recip2 <- 
  recipe(y ~ ., data = df) |> 
  step_nearmiss(y, under_ratio = 1.1) |> 
  step_smote(y, over_ratio = 0.5) |> 
  step_normalize(all_predictors())

df |> ggplot(aes(y)) + geom_bar() + coord_flip()
prep(recip1) |> juice() |> ggplot(aes(y)) + geom_bar() + coord_flip()
prep(recip2) |> juice() |> ggplot(aes(y)) + geom_bar() + coord_flip()

