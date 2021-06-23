library(tidyverse)
library(furrr)
library(tidymodels)

# TODO EDA of prepped data

test <- read_rds('./data/test_df_scored.rds') 
train <- read_rds('./data/train_df_scored.rds')

sum(is.na(test))

sum(is.na(train))

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
  count(subfamily) |> 
  print(n=20)


classed <- 
  test |> 
  select(acc, subfamily, Arch1:Xer) |> 
  pivot_longer(cols = Arch1:Xer, names_to = 'hmm_name', values_to = 'hmm_score') |> 
  
  group_by(acc, subfamily) |> 
  filter(hmm_score == max(hmm_score)) |> 
  left_join(gathering_thresholds) |> 
  mutate(hmm_name = ifelse(hmm_score < threshold, 'non_integrase', hmm_name)) |> 
  ungroup()

  
classed |> 
  count(subfamily, hmm_name) |> 
  left_join(
    test |> count(subfamily, name = 'count')
    ) |> 
  mutate(predicted = hmm_name, 
         actual = fct_rev(subfamily),
         pct_of_actual = n*100 / count) |> 
  ggplot(aes(predicted, actual, fill = pct_of_actual, color  = pct_of_actual)) +
  geom_point(aes(size = n)) +
  scale_size(trans = 'log10') +
  scale_fill_viridis_c() +
  scale_color_viridis_c() +
  theme_light() +
  theme(axis.text.x.bottom = element_text(angle = -90,hjust = 0, vjust = 0))  

  
take_best_hmm <- function(df){
  
}
