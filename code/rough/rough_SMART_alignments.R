# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("Biostrings")
# BiocManager::install("DECIPHER")

library(Biostrings)
library(DECIPHER)
library(tidyverse)

# make dataframe of fastas for each subfamily
subfamilies <- 
  tibble(path = Sys.glob('./data/SMART/domain_fasta/*.fasta')) |> 
  # read, parse fasta
  mutate(name = str_remove_all(path, './data/SMART/domain_fasta/|\\.fasta'),
         fasta = map(path, readAAStringSet),
         nseq = map_int(fasta, ~length(.x)),
         prot_names = map(fasta, ~names(.x)),
         prot_seqs = map(fasta, paste)
         ) |> 
  # join SMART ids
  left_join(read_csv('./data/SMART_YR_ids.csv') |> janitor::clean_names())  |>  
  select(name, id, nseq, everything()) |> 
  # put smallest families first
  arrange(nseq)

glimpse(subfamilies)


# Create alignments using decipher from fasta file, given path
make_alignment <- function(name, path){
  message(paste("\n\nAligning subfamily:", name))
  aligned <- AlignSeqs(readAAStringSet(path))
  outpath <- paste0('./data/SMART/domain_alignments/', name, '.aln')
  writeXStringSet(aligned, filepath = outpath)
  return(aligned)
}

# Do all the alignments. Will create aligned fasta (.aln) files in 
# './data/SMART/domain_alignments/'
# Takes many hours!!
subfamilies <- subfamilies |> 
  mutate(aligned = map2(.x = name, .y = path, .f = ~make_alignment(.x,.y)))

