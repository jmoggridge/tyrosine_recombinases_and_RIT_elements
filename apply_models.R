## apply classifiers to protein sequence...
library(tidyverse)
library(progress)
library(glue)
library(Biostrings)
library(furrr)
library(tidymodels)

# classifiers
knn_model <- read_rds('./models/knn_classifier.rds')
glmnet_model <- read_rds('./models/glmnet_classifier.rds')

# some data that is necessary for classify_proteins()'s hmmsearches
hmmsearch_filler <- read_rds('./data/hmmsearch_filler.rds')

read_fasta_to_tibble <- function(fasta_file){
  fa <- Biostrings::readAAStringSet(fasta_file)
  tibble(name = names(fa),
         protein = paste(fa))
}
# read_fasta_to_tibble('data/test_example.fa')

# classify proteins from unnested genbank feature table,
# returns the integrase classification from the model
classify_proteins <- function(df, name, protein){
  
  # make temp fasta file of seqs to score against hmm library
  # include filler sequences at head of fasta file
  fasta_path <- tempfile()
  proteins <- df |> pull({{protein}})
  sequences <- c(hmmsearch_filler$seq, ft$translation)
  fasta <- Biostrings::AAStringSet(sequences)
  names(fasta) <- c(hmmsearch_filler$name, ft$protein_id)
  Biostrings::writeXStringSet(fasta, fasta_path)  
  
  # setup paths for hmms and table outputs; 
  # build hmmsearch calls (but not executed until next step)
  temp_dir <- tempdir()
  junk <- tempfile()
  hmmer_call <- glue('hmmsearch --noali -o {junk} --tblout')
  hmmsearches <- 
    tibble(hmm_path = Sys.glob(glue('./models/HMM/*'))) |> 
    mutate(
      out_path = hmm_path |> 
        str_replace(glue('./models/HMM/'), glue('{temp_dir}/')) |> 
        str_replace('\\.hmm', '.tbl'),
      calls = glue('{hmmer_call} {out_path} {hmm_path} {fasta_path}')
    )
  
  # do hmmsearch calls in parallel
  plan(multisession, workers = availableCores())
  hmmsearches <- hmmsearches |> 
    mutate(hmmsearches = future_map_dbl(calls, ~ system(.x))) 
  
  # collect hmmsearch results from temp_dir
  ft_scored <- ft |> 
    mutate(acc = protein_id) |> 
    join_hmmsearches2(files = hmmsearches$out_path) |>
    relocate(-c(Arch1:Xer)) |> 
    arrange(feat_id)
  
  # make integrase prediction with classifiers, find consensus
  ft_classed <- ft_scored |> 
    bind_cols(
      knn_pred = predict(knn_model, new_data = ft_scored) |> 
        pull(.pred_class),
      glmnet_pred = predict(glmnet_model, new_data = ft_scored) |> 
        pull(.pred_class)
    ) |> 
    select(-c(Arch1:Xer)) |> 
    mutate(consensus_pred = case_when(
      knn_pred == glmnet_pred ~ paste(knn_pred),
      TRUE ~ 'no consensus'
    ))
}
