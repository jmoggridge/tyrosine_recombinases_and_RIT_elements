# Part 3 - Step 3: Classify the proteins from phage genomes

library(tidyverse)
library(glue)
source('R/P3_predict_functions.R')

phage <- read_rds('data/ncbi_viral_genomes/genbank_parsed.rds')
glimpse(phage)

phage_meta <- read_rds('data/ncbi_viral_genomes/phage_meta.rds')

# unnest genes & their translations; make an id for each
phage_proteins <- phage |>
  select(accession, feature_table) |> 
  unnest(feature_table) |> 
  filter(!is.na(translation)) |> 
  mutate(prot_id = glue::glue('{accession}_{feat_id}')) |> 
  select(prot_id, translation)


classifier_rs <- 
  phage_proteins |>
  predict_integrases(names = prot_id, seqs = translation) |> 
  # should take this out; some names are messed up 
  select(-contains('.y'), -contains('.x')) |> 
  dplyr::rename(accession = acc, prot_id = name, translation = seq)

write_rds(classifier_rs, 'data/ncbi_viral_genomes/classifier_results.rds')
beepr::beep()
  
