### Multiple sequence alignment of subfamilies

# I aligned the domain AA sequences for each subfamily with the decipher package (with default settings). This creates the aligned fasta files in `./data/SMART/domain_alignments` as they are finished. The Int_SXT, Int_Tn916, and Xer sub-families take many hours to align.
# BiocManager::install("DECIPHER")
library(tidyverse)
library(DECIPHER)

# Create alignment of for each subfamily using decipher; 
# reads fasta files from glob of paths
# writes alignments as they finish
make_alignment <- function(name, path){
  message(paste("\n\nAligning subfamily:", name))
  aligned <- AlignSeqs(readAAStringSet(path))
  # outpath <- paste0('./data/SMART/domain_alignments/', name, '.aln')
  # writeXStringSet(aligned, filepath = outpath)
  return(aligned)
}

# Do all the alignments. 
# Creates .aln files in ./data/SMART/domain_alignments 
# Don't need to save obj when finished
# Takes 1 day or more
domains <- read_rds('./data/SMART_db.rds')
domains <- domains |> 
  mutate(aligned = map2(.x = subfamily, 
                        .y = dom_path,
                        .f = ~make_alignment(.x,.y)))
