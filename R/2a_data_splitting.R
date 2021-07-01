## Data splitting for classfier

# Performs a nested CV where the training data are used to create domain subfamily alignments and train HMMs. Then the HMMs aer used to score train and test sequences. The scores are used to predict the class of the test data, and performance metrics for the different models are compared in the inner CV (model tuning) and outer CV (model selection).


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
  # TODO decide on downsample...
  # 10k per subfamily
  group_by(subfamily) |> 
  # slice_sample(n = 10000, replace = F) |> 
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

test <- testing(df_split)
test |> count(subfamily) |> print.AsIs()

write_rds(train, './data/classif_train_set.rds')
write_rds(test, './data/classif_test_set.rds')

rm(smart_df, non_integrases, non_int_leftout, non_int_downsampled, smart_leftout, smart_downsampled, combined_dataset, leftout_dataset, df_split, train, test)
