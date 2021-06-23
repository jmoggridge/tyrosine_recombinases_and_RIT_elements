library(tidyverse)
library(furrr)
library(tidymodels)

# TODO EDA of prepped data

test <- read_rds('./data/test_df_scored.rds')
train <- read_rds('./data/train_df_scored.rds')

gathering_thresholds <- 
  train |> 
  select(subfamily, Arch1:Xer) |> 
  pivot_longer(cols = Arch1:Xer, names_to = 'hmm_name', values_to = 'hmm_score') |> 
  filter(subfamily == hmm_name) |> 
  group_by(subfamily) |> 
  filter(hmm_score == min(hmm_score)) |> 
  ungroup() |> 
  select(-subfamily)



print.data.frame(gathering_thresholds)
  


test |> 
  select(subfamily, Arch1:Xer) |> 
  pivot
take_best_hmm <- function(df){
  
}
