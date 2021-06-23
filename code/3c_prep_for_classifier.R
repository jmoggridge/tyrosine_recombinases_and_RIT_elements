## 3c_prep_for_classifier.R
# Tidies hmmsearch results for train and test dataframes.
# Parses all the `hmmsearch` results tables and joins to test and train observations.

# Found a couple issues in hmmsearch parsing:
#   - table is space-delimited
#   - need to trim whitespace
#   - the hmm names for assessment hmms have '.train' to remove to make the classifier input columns
#   - need to keep an eye that text isn't being truncated ()
# TODO change 3b to make sure no text is truncated.... (leads to ambiguous acc#)
#   - single acc can have multiple rows; when pivoted longer, have 2 rows to join to train or test df. Only want a single row per protein though, so only want to keep the top score.... Solution: add function `max` to pivot longer to keep only max domain score value & then no lists are created.
# TODO Decide whether to use clean_names for classifier

library(tidyverse)
library(purrr)

## fx ----

# BUG prot_name column is truncated - making some acc numbers not match + get duplicated...

read_hmmersearch_tbl <- function(path){
  # read hmmsearch tblout result
  read_delim(
    file = path,
    na = '-',
    delim = ' ',
    comment = '#', 
    trim_ws = T,
    col_types = cols(),
    col_names =  c(
      "prot_name", "prot_accession", 
      "hmm_name", "hmm_accesion",
      "seq_eval", "seq_score", "seq_bias",
      "best_dom_eval", "best_dom_score", "best_dom_bias",
      'domain_n_exp', 'domain_n_reg', 'domain_n_clu', 'domain_n_ov',
      'domain_n_env', 'domain_n_dom', 'domain_n_rep', 'domain_n_inc',
      "description")
  ) 
}

tidy_hmmsearch <- function(df){
  # keep only the best_domain_score column, 
  # remove '.train' from hmm_name if present
  # return (acc, hmm_name) table with best domain scores as values
  df |> 
    select(prot_name, hmm_name, best_dom_score) |> 
    mutate(hmm_name = str_remove(hmm_name, '.train')) |> 
    distinct() |> 
    pivot_wider(names_from = hmm_name, values_from = best_dom_score) |> 
    unnest(cols = c())
}

parse_hmmsearches <- function(glob){
  # reads all hmmsearch tables in glob
  Sys.glob(glob) |> 
    map(read_hmmersearch_tbl) |> 
    # extracts best domain scores for each subfamily and rectangles data
    map(tidy_hmmsearch) |> 
    reduce(full_join, by = 'prot_name') |> 
    mutate(acc = prot_name) |> 
    select(-prot_name) |> 
    relocate(acc)
}


# # join hmmsearch results to train and test datasets
train <- 
  read_rds('./data/train_df.rds')|>
  left_join(parse_hmmsearches("./data/hmmsearch_res/*.train.tbl"), by = 'acc')

glimpse(train)

test <-
  read_rds('./data/test_df.rds') |> 
  left_join(parse_hmmsearches("./data/hmmsearch_res/*.test.tbl"), by = 'acc')

glimpse(test)


# TODO Add kmer features?

# train <- read_rds('./data/train_df.rds')
# test <- read_rds('./data/test_df.rds')
# tr_search <- parse_hmmsearches("./data/hmmsearch_res/*.train.tbl")
# tr_search |> filter(acc %in% train$acc) 
# tr_search |> filter(!acc %in% train$acc) 
# tr_search |> filter(acc %in% test$acc) 
# tr_search |> filter(!acc %in% test$acc) 
# 
# tst_search<- parse_hmmsearches("./data/hmmsearch_res/*.test.tbl")
# tst_search |> filter(acc %in% test$acc) 
# tst_search |> filter(!acc %in% test$acc) 
# tst_search |> filter(acc %in% train$acc) 
# tst_search |> filter(!acc %in% train$acc) 
# 
