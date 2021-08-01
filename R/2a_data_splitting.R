## Data splitting for classfier

# First, smart data is downsampled to a maximum of 10k per subfamily.
# The non-integrase data is already downsampled to 2k per group (by proteome, name, or domain family).
# Modelling data are joined from Smart and Non-integrase datasets.
# Splits data 75/25 to training and final validation sets for modelling.


## Libraries -----------------------------------------------------------------

library(tidyverse)
library(tidymodels)
library(glue)
library(here)

# MAIN --------------------------------------------------------------------

set.seed(123)

## Combine datasets  ------

# integrase dataset
smart_df <- 
  read_rds('./data/SMART/smart_df.rds') |>
  select(subfamily, acc, prot_seq, dom_seq) |> 
  mutate(subfamily = as_factor(subfamily)) 

smart_downsampled <- smart_df |> 
  # downsample to max of 10k per subfamily
  group_by(subfamily) |> 
  slice_sample(n = 10000, replace = F) |>
  ungroup()

smart_leftout <- anti_join(smart_df, smart_downsampled)

# negative examples (non-integrases)
non_integrases <- 
  read_rds('./data/non_integrase_seqs/nonint_df.rds') |> 
  mutate(dom_seq = NA,
         subfamily = as_factor(subfamily))

non_int_downsampled <- 
  non_integrases |> 
  group_by(group) |> 
  sample_frac(size = 0.5, replace = F) |> 
  ungroup()

non_int_leftout <- 
  anti_join(non_integrases, non_int_downsampled) |> 
  select(subfamily, acc, prot_seq, dom_seq)

non_int_downsampled <- non_int_downsampled |> 
  select(subfamily, acc, prot_seq, dom_seq)


# combine datasets 
combined_dataset <- bind_rows(smart_downsampled, non_int_downsampled) 
leftout_dataset <- bind_rows(smart_leftout, non_int_leftout)

write_rds(combined_dataset, './data/classif_combined_dataset.rds')
write_rds(leftout_dataset, './data/classif_leftout_dataset.rds')

## Initial split -------------------------------------------------

df_split <- initial_split(combined_dataset, 0.75, strata = subfamily)

train <- training(df_split)
train |> count(subfamily) |> print.AsIs()
train |> ggplot(aes(subfamily)) + geom_bar() + coord_flip() + labs(title = 'training data')

test <- testing(df_split)
test |> count(subfamily) |> print.AsIs()
test |> ggplot(aes(subfamily)) + geom_bar() + coord_flip()  + labs(title = 'testing data')


write_rds(train, './data/classif_train_set.rds')
write_rds(test, './data/classif_test_set.rds')


