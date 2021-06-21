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
plan(multisession, workers = 8)

# load training domains, remove any duplicated sequences
training_domains <- read_rds('./data/SMART/smart_train.rds') |> 
  group_by(dom_seq) |> 
  sample_n(1) |> 
  ungroup() |> 
  select(subfamily, acc, dom_seq)
  
# subset to check code works
training_domains <- training_domains |> 
  group_by(subfamily) |> 
  slice_head(n=10) |> 
  ungroup()

do_alignment <- function(df, dest){
  print(paste("Working on subfamily:",  df$subfam[1]))
  outpath <- paste0(dest, df$subfam[1], '.aln')
  aa_set <- Biostrings::AAStringSet(df |> pull(dom_seq))
  names(aa_set) <- df$acc
  aligned <- AlignSeqs(aa_set)
  writeXStringSet(aligned, filepath = outpath)
  return(aligned)
}
aligns <- training_domains |> 
  mutate(subfam = subfamily) |> 
  group_by(subfamily) |> 
  nest() |> 
  mutate(aligned = future_map(data, ~do_alignment(.x, dest = './data/SMART/')))

aligns |> 
  select(subfamily, aligned) |> 
  mutate(acc = map(aligned, names),
         aln = map(aligned, paste)) |> 
  unnest(cols = c(acc, aln)) |> 
  View()


# # function to perform alignment and save a fasta file
# do_alignment <- function(df, dest = './data/SMART/'){
#   print(paste("Working on subfamily:",  df$subfam[1]))
#   outpath <- paste0(dest, df$subfam[1], '.aln')
#   aa_set <- Biostrings::AAStringSet(df |> pull(dom_seq))
#   names(aa_set) <- df$acc  
#   aligned <- AlignSeqs(aa_set)
#   writeXStringSet(aligned, filepath = outpath)
#   return(aligned)
# }

# function to perform alignment and save a fasta file
do_alignment <- function(df, dest = './data/SMART/'){
  print(paste("Working on subfamily:",  df$subfam[1]))
  outpath <- paste0(dest, df$subfam[1], '.aln')
  aa_set <- Biostrings::AAStringSet(df |> pull(dom_seq))
  names(aa_set) <- df$acc  
  aligned <- AlignSeqs(aa_set)
  writeXStringSet(aligned, filepath = outpath)
  return(aligned)
}



library(furrr)

training_domains |> 
  # need to duplicate group labels for filenames
  mutate(subfam = subfamily) |> 
  group_by(subfamily) |> 
  group_walk(.f = ~do_alignment(.x, dest = './data/SMART/domain_align_training/'))

# # create set fasta files of training domains for alignment
# system('mkdir ./data/SMART/training_domain_fasta')
# smart_train |> 
#   mutate(label = subfamily) |> 
#   group_by(subfamily) |> 
#   group_walk(.f = ~make_fasta(.x))
# 
# 
#   
# # Create alignment of for each subfamily using decipher; 
# # writes alignments as they finish
# make_alignment <- function(name, path){
#   message(paste("\n\nAligning subfamily:", name))
#   aligned <- AlignSeqs(readAAStringSet(path))
#   outpath <- paste0('./data/SMART/domain_alignments/', name, '.aln')
#   writeXStringSet(aligned, filepath = outpath)
#   return(aligned)
# }
# 
