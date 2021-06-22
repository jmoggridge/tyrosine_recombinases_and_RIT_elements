## 3c_prep_for_classifier.R
# Tidies hmmsearch results for train and test dataframes.
# Parses all the `hmmsearch` results tables and joins to test and train observations.

library(tidyverse)
library(purrr)

## fx ----
read_hmmersearch_tbl <- function(path){
  # read hmmsearch tblout result
  # keep only the best_domain_score column, 
  # rename it as the subfamily name for the hmm
  read_table(
    file = path,
    na = '-',
    comment = '#',
    col_types = cols(),
    col_names =  c(
      "prot_name", "prot_accession", 
      "hmm_name", "hmm_accesion",
      "seq_eval", "seq_score", "seq_bias",
      "best_dom_eval", "best_dom_score", "best_dom_bias",
      'domain_n_exp', 'domain_n_reg', 'domain_n_clu', 'domain_n_ov',
      'domain_n_env', 'domain_n_dom', 'domain_n_rep', 'domain_n_inc',
      "description")
  ) |> 
    select(prot_name, hmm_name, best_dom_score) |> 
    distinct() |> 
    pivot_wider(names_from = hmm_name, values_from = best_dom_score) |> 
    unnest(cols = c())
}

# collects all hmmsearch tables from glob, join cols by acc#, and return df.
parse_hmmsearches <- function(glob){
  Sys.glob(glob) |> 
    map(read_hmmersearch_tbl) |> 
    purrr::reduce(full_join, by = 'prot_name') |> 
    mutate(acc = prot_name) |> 
    select(-prot_name)
}

# TODO CONTINUE THIS SCRIPT


# join hmmsearch results to train and test datasets
train <- 
  read_rds('./data/train_df.rds')
  left_join(parse_hmmsearches("./data/hmmsearch_res/*.train.tbl"), by = 'acc')
  
test <- 
  read_rds('./data/test_df.rds')
  left_join(parse_hmmsearches("./data/hmmsearch_res/*.test.tbl"), by = 'acc')






