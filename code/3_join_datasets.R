### Join integrase and non-integrase data ####

library(tidyverse)
library(rsample)
library(Biostrings)


# # non_integrasess
# non_integrases <- read_rds('./data/non_integrase_seqs/non_integrases_df.rds') 
# glimpse(non_integrases)


# smart data
smart_train <- read_rds('./data/SMART/smart_train.rds')
smart_test <- read_rds('./data/SMART/smart_test.rds')
glimpse(smart_test)
glimpse(smart_train)

# join training datasets
train_df <- 
  bind_rows(
    smart_train |> select(subfamily, acc, description, prot_seq),
    nonint_train |> select(subfamily, acc, description, prot_seq)
  )
glimpse(train_df)

# join testing datasets
test_df <- 
  bind_rows(
    smart_test |> select(subfamily, acc, description, prot_seq),
    nonint_test |> select(subfamily, acc, description, prot_seq)
  )
glimpse(test_df)

# save dataframes for joining hmmsearch scores with
write_rds(train_df, './data/train_df.rds', compress = 'gz')
write_rds(test_df, './data/test_df.rds', compress = 'gz')

## create fasta file for hmmsearch scoring vs 20 domain HMMs
write_fasta <- function(df, dest){
  fasta <- Biostrings::AAStringSet(df$prot_seq)
  names(fasta) <- df$acc
  Biostrings::writeXStringSet(fasta, filepath = dest)
}

write_fasta(train_df, './data/train_seq.fa')
write_fasta(test_df, './data/test_seq.fa')


## Proceed to scoring seqs with hmmsearch: 