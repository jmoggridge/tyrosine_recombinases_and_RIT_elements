## 3c_prep_for_classifier.R
# Tidies hmmsearch results for train and test dataframes.
# Parses all the `hmmsearch` results tables and joins to test and train observations.

# Found a couple issues in hmmsearch parsing:
#   - table is space-delimited
#   - need to trim whitespace
#   - the hmm names for assessment hmms have '.train' to remove to make the classifier input columns
#   - need to keep an eye that text isn't being truncated ()
#   - single acc can have multiple rows; when pivoted longer, have 2 rows to join to train or test df. Only want a single row per protein though, so only want to keep the top score.... Solution: add function `max` to pivot longer to keep only max domain score value & then no lists are created.

# TODO Decide whether to use clean_names for classifier

library(tidyverse)
library(purrr)
library(Biostrings)
library(kmer)


## ---- fx ----

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
      "prot_name", "prot_accession", "hmm_name", "hmm_accesion", "seq_eval",
      "seq_score", "seq_bias", "best_dom_eval", "best_dom_score", "best_dom_bias",
      'domain_n_exp', 'domain_n_reg', 'domain_n_clu', 'domain_n_ov', 'domain_n_env',
      'domain_n_dom', 'domain_n_rep', 'domain_n_inc', "description")
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
    purrr::reduce(full_join, by = 'prot_name') |> 
    mutate(acc = prot_name) |> 
    select(-prot_name) |> 
    relocate(acc)
}

as_AAbin <- function(z, simplify = FALSE){
# need to convert sequence from character to AAbin object
  res <- if(length(z) == 1 & simplify) charToRaw(z) else lapply(z, charToRaw)
  attr(res, "rerep.names") <- attr(z, "rerep.names")
  attr(res, "rerep.pointers") <- attr(z, "rerep.pointers")
  class(res) <- "AAbin"
  return(res)
}

add_2mers <- function(df, seq){
# converts seq to AAbin obj, computes kmer proportions
  kmers <- df |> 
    pull({{seq}}) |> 
    as_AAbin() |> 
    # count kmers; take proportion by normalizing counts across each row
    kcount(k = 2, residues = 'AA') |> 
    apply(1, function(x) x/sum(x)) |> 
    # matrix to df
    t() |> 
    as_tibble()
  df |> bind_cols(kmers)
}

## IDEA add dayhoff 5mers

# add_5mers <- function(df, seq){
#   # as for add_2mers but with 5mers over reduced dayhoff alphabet
#   kmers <- df |> 
#     pull({{seq}}) |> 
#     as_AAbin() |> 
#     # count kmers; take proportion by normalizing counts across each row
#     kcount(k = 5, residues = 'AA') |> 
#     apply(1, function(x) x/sum(x)) |> 
#     # matrix to df
#     t() |> 
#     as_tibble()
#   df |> bind_cols(kmers)
# }

## ---- main ----

# training data
train <- 
  read_rds('./data/train_df.rds') |>
  # join hmmsearch results to dataset
  left_join(
    parse_hmmsearches("./data/hmmsearch_res/*.train.tbl"), 
    by = 'acc'
    ) |> 
  # replace any NAs with zeros
  mutate(across(Arch1:Xer, ~replace_na(.x, 0))) |> 
  # join 2-mer proportions columns
  add_2mers(seq = prot_seq) 

write_rds(train, './data/train_df_scored.rds')
rm(train)

## same thing for test data

test <- 
  read_rds('./data/test_df.rds') |> 
  left_join(
    parse_hmmsearches("./data/hmmsearch_res/*.test.tbl"),
    by = 'acc'
    ) |>
  mutate(across(Arch1:Xer, ~replace_na(.x, 0))) |> 
  add_2mers(seq = prot_seq)

write_rds(test, './data/test_df_scored.rds')
rm(test)



## ----- tests -----

## checks to make sure that all proteins are in the right spot
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


# test |> 
#   select(subfamily, acc, prot_seq) |> 
#   filter(str_detect(prot_seq, '[^ACDEFGHIKLMNPQRSTVWY]')) |> 
#   mutate(wtf_char = str_extract_all(prot_seq, '[^ACDEFGHIKLMNPQRSTVWY]'),
#          wtf_char = map_chr(wtf_char, ~paste0(.x, collapse = ''))) |> 
#   View()
#   
# 
# train |> 
#   select(subfamily, acc, prot_seq) |> 
#   filter(str_detect(prot_seq, '[^ACDEFGHIKLMNPQRSTVWY]')) |> 
#   mutate(wtf_char = str_extract_all(prot_seq, '[^ACDEFGHIKLMNPQRSTVWY]'),
#          wtf_char = map_chr(wtf_char, ~paste0(.x, collapse = ''))) |> 
#   arrange(desc(nchar(wtf_char))) |> 
#   View()
