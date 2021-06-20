library(tidyverse)
library(tidymodels)
library(workflowsets)

set.seed(1)

# make a small dataset with stratified sampling for code development
df <- 
  read_rds('./data/smart_refseqs_hmm_scores.rds') |>
  select(-id) |> 
  relocate(subfamily, contains('prot_')) |> 
  group_by(subfamily) |> 
  sample_frac(0.1) |> 
  ungroup() |> 
  mutate(across(Arch1:Xer, ~replace_na(.x, 0))) |> 
  janitor::clean_names()

# count kmers


glimpse(df)
summary(df)
summary(df$subfamily)


library(corrr)

cor_df <- df |> 
  select(arch1:xer) |> 
  mutate(across(everything(), ~(.x - mean(.x)) / sd(.x))) |> 
  correlate() |> 
  pivot_longer(-term, names_to = 'term2', values_to = 'cor') |> 
  mutate(cor = ifelse(is.na(cor), 0, cor)) |> 
  arrange(-cor) |> 
  print(n=50)

cor_df |> 
  group_by(term) |> 
  summarize(mean = mean(cor)) |> 
  arrange(-mean)

