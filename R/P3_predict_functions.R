
require(tidyverse)
require(tidymodels)
require(glue)
require(furrr)
require(Biostrings)

# load models
models <- 
  Sys.glob('models/*_classifier.rds') |> 
  purrr::map(read_rds)
names(models) <- c('glmnet', 'knn', 'rf')

# path to hmms
hmm_folder <- 'models/hmm/'

# TODO make the filler unnecessary by checking if hmmsearch res 
# file is empty with ifelse statement before trying to parse it...

# filler data that is necessary for classify_proteins()'s hmmsearches
# such that they don't return an empty file...
hmmsearch_filler <- readr::read_rds('models/hmmsearch_filler.rds')

## functions to read & join hmmsearch output for 20 HMMs
read_hmmsearch <- function(path){
  # read hmmsearch tblout result for one dataset vs one HMM; 
  # return df of (acc, <hmm_name>)
  suppressWarnings(
    hmmsearch_rs <- 
      read_delim(
        file = path, na = '-', delim = ' ', comment = '#',  trim_ws = T,
        col_types = cols(), col_names = FALSE, 
      ) |> 
      transmute(acc = X1, hmm_name = X3, best_dom_score = X9) |> 
      # keep only the best hmm score for each protein
      group_by(acc) |> filter(best_dom_score == max(best_dom_score)) |> 
      ungroup() |>  distinct() |> 
      # name the values column after the subfamily hmm
      pivot_wider(names_from = hmm_name, values_from = best_dom_score) 
  )
  
  return(hmmsearch_rs)
}

join_hmmsearches <- function(df, files){
  # df must have identifiers column to match hmmer output 
  # (need to match the headers in the fasta passed to hmmsearch)
  # read and combine hmmsearch data with reduce
  # join scores to input data & replace any NAs with zeros
  files |> 
    map(read_hmmsearch) |> 
    purrr::reduce(.f = full_join, by = 'acc') |> 
    right_join(df, by = 'acc') |>
    mutate(across(Arch1:Xer, ~replace_na(.x, 0))) |> 
    relocate(everything(), Arch1:Xer)
}


do_hmmsearches <- function(df, names, seqs){}


# classify proteins from unnested genbank feature table,
# returns the integrase classification from the model
predict_integrases <- function(df, names, seqs){
  
  # get sequences and their names into vectors to create one fasta
  seq_v <- df |> pull({{seqs}})
  name_v <- df |> pull({{names}})
  

  # df needs to have an id column called 'acc' for joining hmmsearch results
  # to work models to work
  df <- df |> 
    mutate(acc = word(name_v, 1L),
           name = {{names}},
           seq = {{seqs}}) |> 
    select(-{{names}}, -{{seqs}})
  
  # make temp fasta file of seqs to score against hmm library
  # include filler sequences first
  fasta_path <- tempfile()
  fasta <- Biostrings::AAStringSet(c(hmmsearch_filler$seq, seq_v))
  names(fasta) <- c(hmmsearch_filler$name, name_v)
  Biostrings::writeXStringSet(fasta, fasta_path)  
  
  # setup paths for hmms and table outputs; 
  # build hmmsearch calls (but not executed until next step)
  temp_dir <- tempdir()
  junk <- tempfile()
  plan(multisession, workers = availableCores())
  hmmer_call <- glue::glue('hmmsearch --noali -o {junk} --tblout')
  
  # do hmmsearch calls in parallel
  hmmsearches <- 
    tibble(hmm_path = Sys.glob(glue::glue('models/hmm/*'))) |> 
    mutate(
      out_path = hmm_path |> 
        str_replace(
          glue::glue('models/hmm/'), 
          glue::glue('{temp_dir}/')) |> 
        str_replace('\\.hmm', '.tbl'),
      calls = glue::glue('{hmmer_call} {out_path} {hmm_path} {fasta_path}')
    ) |> 
    mutate(hmmsearches = future_map_dbl(calls, ~ system(.x))) 
  
  
  # collect hmmsearch results from temp_dir and join to dataframe
  df_scored <- 
    df |> 
    join_hmmsearches(files = hmmsearches$out_path) |>
    relocate(-c(Arch1:Xer))
  
  # make integrase prediction with classifiers, find consensus
  df_classed <-
    df_scored |> 
    bind_cols(
      glmnet_pred = predict(models$glmnet, new_data = df_scored) |> 
        pull(.pred_class),
      knn_pred = predict(models$knn, new_data = df_scored) |> 
        pull(.pred_class),
      rf_pred = predict(models$rf, new_data = df_scored) |> 
        pull(.pred_class),
    ) |> 
    mutate(consensus_pred = case_when(
      knn_pred == glmnet_pred ~ paste(knn_pred),
      knn_pred == rf_pred ~ paste(knn_pred),
      glmnet_pred == rf_pred ~ paste(glmnet_pred),
      TRUE ~ 'no consensus'
    ))
  
  return(
    df |> left_join(df_classed, by = c("name", "seq", "acc"))
    )
}


