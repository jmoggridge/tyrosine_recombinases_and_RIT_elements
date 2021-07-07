library(tidyverse)
library(Biostrings)
library(tidymodels)
library(furrr)
library(tictoc)

# TODO FINAL CLASSIFIER FITTING STEP


combined_dataset <- 
  read_rds('./data/classif_combined_dataset.rds')
combined_dataset
combined_dataset |> count(subfamily) |> print(n=25)
skimr::skim(combined_dataset)

# TODO full align
# aligns (use code from earlier)

# TODO full hmmm

# TODO full scores

# score a dataframe of sequences against the HMM library with hmmsearch
# takes protein seq column and makes temporary fasta file to send to hmmer
hmmsearch_scores2 <- function(df, hmm_path, out_folder, tag){
  
  # make temp fasta file of seqs to score against hmm library
  fasta_path <- tempfile()
  fasta <- Biostrings::AAStringSet(df$prot_seq)
  names(fasta) <- df$acc
  writeXStringSet(fasta, fasta_path)  
  
  junk <- tempfile()
  # setup paths for hmms and table outputs; build hmmsearch calls
  hmmsearches <- 
    tibble(hmm_path = Sys.glob(glue('{hmm_path}*'))) |> 
    mutate(
      out_path = hmm_path |> 
        str_replace('/hmm/', '/hmmsearch/') |> 
        str_replace('\\.hmm', glue('.tbl')),
      calls = glue('hmmsearch --noali -o {junk} --tblout {out_path} {hmm_path} {fasta_path}')
    )
  # do hmmsearch calls in parallel
  plan(multisession, workers = availableCores())
  hmmsearches <- hmmsearches |> 
    mutate(hmmsearches = future_map_dbl(calls, ~ system(.x))) 
  return(hmmsearches)
}

hmm_scores <- 
  hmmsearch_scores2(
    df = combined_dataset, 
    hmm_path = './data/SMART/domain_hmm_final/',
    )

  


# TODO FIT IT

# fit final model: 10-NN trained on combined dataset
# specify modelling workflow with scaling and ignore seqs and id
recip <- 
  recipe(subfamily ~ ., data = combined_dataset) |> 
  update_role(c('acc', contains('seq')), new_role = "id") |> 
  step_scale(all_predictors())
knn_spec <- 
  nearest_neighbor(mode = 'classification', neighbors = 10) |> 
  set_engine('kknn') 
knn_wkfl <- workflow() |> 
  add_recipe(recip) |> 
  add_model(knn_spec) |> 
  fit(data = combined_dataset)
  
  

# given df with columns acc & prot_seq, and path to 20 subfamily HMMs
# score seqs against hmms -save to temp folder,
# parse hmmer files and return df with HMM scores
score_seqs <- function(df, acc, prot_seq, hmm_folder){
  
}
  
  