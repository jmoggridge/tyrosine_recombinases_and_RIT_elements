
## join integrase and non-integrase data

library(tidyverse)
library(rsample)

# smart data
smart_train <- read_rds('./data/SMART/smart_train.rds')
smart_test <- read_rds('./data/SMART/smart_test.rds')

# non_integrases
non_integrase <- read_rds('./data/non_integrases_df.rds') |> 
  select(-title) |> 
  rename(group = name,
         prot_seq = seq,
         acc = prot_name,
         description = prot_description) |> 
  mutate(subfamily = 'non_integrase')

glimpse(non_integrase)

non_integrase |> 
  pull(prot_seq) |> 
  toupper() |> 
  str_count('X') |> 
  summary()

nonint_split <- initial_split(non_integrase)
nonint_train <- 


## create fasta file for HMM scoring
# 
# # create fasta files for scoring train and test_data
# dest <- './data/train_seqs.fa'
# fasta <- Biostrings::AAStringSet(seq)
# names(fasta) <- acc
# writeXStringSet(fasta, filepath = dest)

## score HMMs





