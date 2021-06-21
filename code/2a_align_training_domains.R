## Alignments
# one set for training HMMs for classification assessment
# one set for full data HMMs for final classifier

library(tidyverse)
library(furrr)
library(Biostrings)
library(DECIPHER)


# pass a dataframe with subfamily, acc, dom_seq
# return same df with alignment stringset column
do_alignment <- function(df, dest){
  print(paste("Working on subfamily:",  df$subfam[1]))
  outpath <- paste0(dest, df$subfam[1], '.aln')
  aa_set <- Biostrings::AAStringSet(df |> pull(dom_seq))
  names(aa_set) <- df$acc
  aligned <- AlignSeqs(aa_set)
  writeXStringSet(aligned, filepath = outpath)
  return(aligned)
}

### Make training domain alignments ----
system('mkdir ./data/SMART/domain_align_training/')

# load training domains, remove any duplicated sequences
training_domains <- read_rds('./data/SMART/smart_train.rds') |> 
  group_by(dom_seq) |> 
  sample_n(1) |> 
  ungroup() |> 
  select(subfamily, acc, dom_seq)
  
## # subset to check code works
## training_domains <- training_domains |> 
##   group_by(subfamily) |> 
##   slice_head(n=1000) |> 
##   ungroup()

# apply alignment to each subfamily of domains
plan(multisession, workers = 8)
aligns <- training_domains |> 
  mutate(subfam = paste0(subfamily, '.train')) |> 
  group_by(subfamily) |> 
  nest() |> 
  ungroup() |> 
  mutate(aligned = future_map(data, ~do_alignment(.x, dest = './data/SMART/domain_align_training/')))

aligns <- aligns |> 
  select(subfamily, aligned) |> 
  mutate(acc = map(aligned, names),
         dom_aln = map(aligned, paste)) |> 
  unnest(cols = c(acc, aln)) |> 
  glimpse()

training_df <- full_join(training_domains, aligns, by = c("subfamily", "acc")) |> 
  full_join(read_rds('./data/SMART/smart_train.rds'))

write_rds(training_df, './data/smart_train_aligned.rds')

## send alignments to build HMMs...



