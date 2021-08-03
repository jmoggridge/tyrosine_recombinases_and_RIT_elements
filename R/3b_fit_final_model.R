library(tidyverse)
library(Biostrings)
library(tidymodels)
library(furrr)
library(tictoc)

source('./R/00_functions.R')
dir.create('./models/')

combined_dataset <- 
  read_rds('./data/classif_combined_dataset.rds')
combined_dataset
combined_dataset |> count(subfamily) |> print(n=25)
skimr::skim(combined_dataset)

# TODO full align from earlier script

# combined dataset |> align domains

# TODO full hmmms from earlier script
# combined dataset |> build hmms


# score a dataframe of sequences against the HMM library with hmmsearch
# takes protein seq column and makes temporary fasta file to send to hmmer
hmmsearch_scores2 <- function(df, hmm_folder, out_folder){
  
  # make temp fasta file of seqs to score against hmm library
  fasta_path <- tempfile()
  fasta <- Biostrings::AAStringSet(df$prot_seq)
  names(fasta) <- df$acc
  writeXStringSet(fasta, fasta_path)  
  
  junk <- tempfile()
  # setup paths for hmms and table outputs; build hmmsearch calls
  hmmsearches <- 
    tibble(hmm_path = Sys.glob(glue('{hmm_folder}*'))) |> 
    mutate(
      out_path = hmm_path |> 
        str_replace(hmm_folder, out_folder) |> 
        str_replace('\\.hmm', glue('.tbl')),
      calls = glue('hmmsearch --noali -o {junk} --tblout {out_path} {hmm_path} {fasta_path}')
    )
  # do hmmsearch calls in parallel
  plan(multisession, workers = availableCores())
  hmmsearches <- hmmsearches |> 
    mutate(hmmsearches = future_map_dbl(calls, ~ system(.x))) 
  return(hmmsearches)
}

# dir to put the HMM scores
dir.create( './data/SMART/final_hmmscores/')

# Score all sequences in the combined dataset
hmm_scores <- 
  hmmsearch_scores2(
    df = combined_dataset, 
    hmm_folder = './data/SMART/domain_hmm_final/',
    out_folder = './data/SMART/final_hmmscores/'
    )

# scoring data for all sequences, not normalized or smote'd
combined_prepped <- 
  join_hmmsearches(combined_dataset, hmm_scores$out_path) |> 
  select(acc, subfamily, Arch1:Xer)

write_rds(combined_prepped, './data/classif_combined_prepped_dataset.rds', 
          compress = 'gz')
rm(combined_dataset, hmm_scores)


# FIT MODELS -----

# fit final model: 10-NN trained on combined dataset
# specify modelling workflow with scaling and ignore seqs and id
recip <- 
  recipe(subfamily ~ ., data = combined_prepped) |> 
  update_role(acc, new_role = "id") |> 
  step_smote(subfamily, over_ratio = 0.25) |> 
  step_normalize(all_predictors())
recip

# normalized, smote'd HMM scores data for the combined dataset
write_rds(
  x = recip |> prep() |> juice(),
  file = './data/classif_combined_prepped_normalz_smote_dataset.rds',
  compress = 'gz' 
)
best_models <- 
  read_rds('./results/final_model_07-11/final_validation_results.rds') |> 
  left_join(read_rds('./data/unfitted_parsnip_model_set.rds'))
best_models

## Model specifications to finalize ----

glmnet_spec <- best_models |> 
  filter(model_type == 'multinomial regression') |> 
  unnest(params) |> 
  filter(penalty == max(penalty)) |> 
  pull(spec) |> 
  pluck(1)
glmnet_spec

knn_spec <- best_models |> 
  filter(model_type == 'nearest neighbor') |> 
  pull(spec) |> 
  pluck(1)
knn_spec

rf_spec <- best_models |> 
  filter(model_type == 'random forest') |> 
  pull(spec) |> 
  pluck(1) 
rf_spec

rm(best_models)

## fit model workflowss -------

glmnet_classifier <- 
  workflow() |> 
  add_recipe(recip) |> 
  add_model(glmnet_spec) |> 
  fit(data = combined_prepped)
print(glmnet_classifier)
write_rds(glmnet_classifier, './models/glmnet_classifier.rds')
rm(glmnet_classifier)

knn_classifier <- 
  workflow() |> 
  add_recipe(recip) |> 
  add_model(knn_spec) |> 
  fit(data = combined_prepped)
print(knn_classifier)
write_rds(knn_classifier, './models/knn_classifier.rds')
rm(knn_classifier)

rf_classifier <-  
  workflow() |> 
  add_recipe(recip) |> 
  add_model(rf_spec) |> 
  fit(data = combined_prepped)
print(rf_classifier)
write_rds(rf_classifier, './models/rf_classifier.rds')
rm(rf_classifier)


## rule-based model -----

# unnormalized threshold
thresholds <-
  combined_prepped |> 
  select(subfamily, Arch1:Xer) |> 
  pivot_longer(Arch1:Xer, names_to = 'HMM', values_to = 'score') |> 
  filter(HMM == subfamily) |> 
  group_by(HMM) |> 
  filter(score == min(score)) |> 
  transmute(HMM, threshold = score)

thresholds
write_rds(thresholds, './models/hmm_thresholds.rds')

# normalized + smoted thresholds
normalized_and_smoted <- 
  recip |> prep() |> juice()

normalized_and_smoted |> count(subfamily) |> ggplot(aes(n, subfamily)) + geom_col()

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

beepr::beep()



## OLD
# recip <- 
#   recipe(subfamily ~ ., data = combined_prepped) |> 
#   update_role(acc, new_role = "id") |> 
#   step_scale(all_predictors())
