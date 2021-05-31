library(tidyverse)
library(Biostrings)

## Prokaryotic viral orthologous groups database
# 2016 update version
# Viruses that infect bacteria or archaea
# Formely known as Phage Orthologous Groups
# 9518 orthologous groups
# nearly 3k complete genomes
# 296,595 proteins



pvog_raw <- read_delim("./Data/pVOGs_AllFamilyProteinList.tsv", '\t') %>%
  janitor::clean_names()

names(pvog_raw) <- c('vog_number', 'n_proteins', 'protein_accessions')
  
# remove divider lines
# find family headers
# make it so that each row has the name of the virus family
# remove the header rows with family name
# remove hashes from VOG_number
pvog <- pvog_raw %>% 
  filter(!str_detect(vog_number, '#===')) %>% 
  mutate(family = ifelse(
    str_detect(vog_number, '^#'), NA, vog_number)
    ) %>% 
  fill(family) %>% 
  filter(str_detect(vog_number, '^#'))

glimpse(pvog)

# expand protein accession list
pvog %>% 
  mutate(protein_accessions = str_split(protein_accessions, ','))











