library(Biostrings)
library(tidyverse)
# uniprot annotations retrieval for SMART YR domain proteins

Sys.glob("./Data/SMART_tyrosine_recombinase_proteins/*")


RitC_fasta <- readBStringSet("./Data/SMART_tyrosine_recombinase_proteins/RitC_proteins.fasta") 

RitC <- tibble(
  name = names(RitC_fasta),
  seq = paste(RitC_fasta)
)

RitC$name[1]
RitC$seq[1]
