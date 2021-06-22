## Alignments
# one set for training HMMs for classification assessment
# one set for full data HMMs for final classifier

## need to recheck code ~ filenames / paths may have been changed...
library(tidyverse)
library(furrr)
library(Biostrings)
library(DECIPHER)

## Align seqs fx -------
# pass a dataframe with subfamily, acc, dom_seq
# return same df with alignment stringset column
do_alignment <- function(df, dest){
  print(paste('starting...', df$subfam[1]))
  outpath <- paste0(dest, df$subfam[1], '.aln')
  aa_set <- Biostrings::AAStringSet(df |> pull(dom_seq))
  names(aa_set) <- df$acc
  aligned <- AlignSeqs(aa_set, verbose = F, processors = 1)
  writeXStringSet(aligned, filepath = outpath)
  return(aligned)
}


### Organize data by subfamily ----
system('mkdir ./data/SMART/domain_align_allseqs/')

set.seed(123)
all_domains <- 
  # load training domains
  read_rds('./data/SMART/smart_df.rds') |> 
  # remove any duplicated sequences
  group_by(dom_seq) |> 
  sample_n(1) |> 
  ungroup() |> 
  select(subfamily, acc, dom_seq) |>
  # nest names and seqs by family, create subfam label to pass for filename
  mutate(subfam = subfamily) |>
  group_by(subfamily) |>
  nest() |>
  ungroup() |> 
  # arrange by class size to parallelize efficiently
  mutate(nrow = map_int(data, nrow)) |> 
  arrange(desc(nrow)) |> 
  select(-nrow) 

all_domains


## Alignment ---------

# setup parallelization of future_map()
plan(multisession, workers = availableCores() - 1)

# apply alignment to each subfamily of domains
aligns <- all_domains |> 
  mutate(
    aligned = future_map(
      .x = data, 
      .f = ~do_alignment(.x, dest = './data/SMART/domain_align_allseqs/'),
      .options = furrr_options(scheduling = Inf)
    )
  )

aligns
write_rds(aligns, './data/SMART/smart_full_aligned.rds')
