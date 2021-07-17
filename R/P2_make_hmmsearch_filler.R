## P2_make_hmmsearch_filler

# create a dataset to tack onto any queries so that hmmsearch doesn't return a blank file
library(tidyverse)
library(Biostrings)

seqs <- Biostrings::readAAStringSet('./data/hmmsearch_filler_seq.fa')
tibble(name = names(seqs), seq = paste(seqs)) |> 
  write_rds('./data/hmmsearch_filler.rds')
