##

library(tidyverse)
source('./R/00_functions.R')


## Main --------------------------------------------------------------------

knn_model <- read_rds('./results/knn_model_wkfl.rds')
iceberg <- read_tsv('data/iceberg/ICE_db.tsv')

ice <- iceberg |> 
  transmute(acc = prot_accession, prot_seq) |> 
  filter(!is.na(prot_seq)) |> 
  distinct()
ice


# HMM scores

temp_dir <- tempdir()
# make temp fasta file of seqs to score against hmm library
fasta_path <- tempfile()
fasta <- Biostrings::AAStringSet(ice$prot_seq)
names(fasta) <- ice$acc
writeXStringSet(fasta, fasta_path)  

junk <- tempfile()
# setup paths for hmms and table outputs; build hmmsearch calls
hmmsearches <- 
  tibble(hmm_path = Sys.glob(glue('{hmm_path}*'))) |> 
  mutate(
    out_path = hmm_path |> 
      str_replace('./data/SMART/domain_hmm_final/', glue('{temp_dir}/')) |> 
      str_replace('\\.hmm', glue('.tbl')),
    calls = glue('hmmsearch --noali -o {junk} --tblout {out_path} {hmm_path} {fasta_path}')
    )

# do hmmsearch calls in parallel
plan(multisession, workers = availableCores())
hmmsearches <- hmmsearches |> 
  mutate(hmmsearches = future_map_dbl(calls, ~ system(.x))) 

join_hmmsearches2 <- function(df, files){
  # read and combine hmmsearch data
  searches <- map(files, read_hmmsearch) |> 
    purrr::reduce(.f = full_join, by = 'acc') 
  # join scores to input data & replace any NAs with zeros
  df <- df |> 
    left_join(searches, by = 'acc')
  df <- df |>
    mutate(across(Arch1:Xer, ~replace_na(.x, 0))) 
  return(df)
}

ice_scored <- 
  join_hmmsearches2(df = ice, files = hmmsearches$out_path)

preds <- 
  predict(knn_model, new_data = ice_scored) |> 
  bind_cols(ice_scored) |> 
  mutate(prot_accession = acc) |> 
  select(-acc) |> 
  left_join(ice |> mutate(prot_accession = acc))

iceberg_preds <- 
  left_join(iceberg, preds, by = c("prot_seq", "prot_accession"))

summary(iceberg_preds$.pred_class)
iceberg_preds |> 
  count(.pred_class) |> 
  ggplot(aes(n, fct_rev(.pred_class))) +
  geom_col() +
  geom_text(aes(label = n), nudge_x = 800) +
  theme_bw() +
  labs(x = 'count', y = 'YR subfamily prediction', 
       title = 'ICEberg2.0 dataset protein classifications')

iceberg_preds |> 
  filter(str_detect(.pred_class, 'Rit')) |> 
  group_by(type, name, id) |> 
  count(.pred_class) |> 
  ungroup() |> 
  gt::gt() |> 
  gt::tab_options(table.font.size = 11)
