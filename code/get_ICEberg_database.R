library(Biostrings)
library(tidyverse)
 
## The ICEberg2.0 database has a set of 26,566 proteins associated with integrative conjugative genetic elements (ICEs) as those or predicted to be. https://db-mml.sjtu.edu.cn/ICEberg2/download.html. Here, I obtained  data for all experimental-validated elements with intact sequences in DNA and protein fastas for ICEs, IMEs, AICEs, CIMEs, T4SS-type ICEs
#' 

## collect data from ICEberg website
get_ICEberg <- function(element, moltype){
  base <- "https://db-mml.sjtu.edu.cn/ICEberg2/download/"
  exp_intact <- "_experimental_intact.fas"
  link <- paste0(base, element, '_', moltype, exp_intact)
  if (moltype == 'aa'){
    return(Biostrings::readAAStringSet(link))
  } else{
    return(Biostrings::readDNAStringSet(link))
  }
}

# Organize data for all types of elements in single dataframe
MGEs <- 
  tibble(type = c('AICE', 'CIME', 'ICE', 'IME', 'T4SS-type_ICE')) %>% 
  mutate(
    dna_fasta = map(type, ~get_ICEberg(.x, "seq")),
    aa_fasta = map(type, ~get_ICEberg(.x, "aa"))
    )

# save raw nested dataset to file
write_rds(MGEs, "./Data/ICEberg_experimental_intact_database_raw.rds")


# Extract ICE seqs and header labels from fasta (ea. row is an element)
parse_dna_fasta <- function(dna_fasta){
  df <- tibble(
    name = names(dna_fasta),
    dna_seq = paste(dna_fasta)
  ) %>% 
    mutate(header = str_split(name, '\\|'),
           id = map_chr(header, 2),
           name = map_chr(header, 3),
           database = map_chr(header, 4),
           accession = map_chr(header, 5),
           description = map_chr(header, 6),
           start_pos = str_extract(description, "^[0-9]+"),
           end_pos = str_extract(description, "\\.\\.[0-9|]+"),
           end_pos = str_remove_all(end_pos, '\\.'),
    ) %>%
    select(name, id, database, accession, start_pos,
           end_pos, description, dna_seq)
}

dna_df <- MGEs %>% 
  mutate(dna_data = map(dna_fasta, parse_dna_fasta)) %>% 
  select(-contains('fasta')) %>%
  unnest(dna_data)
dna_df


# Extract protein seqs and header labels from fasta (multiple/element)
parse_protein_fasta <- function(aa_fasta){
  df <- tibble(
    header = names(aa_fasta),
    prot_seq = paste(aa_fasta)
  ) %>%
    #   # extract annotation info from sequence headers
    mutate(
      header = map(header, ~unlist(str_split(.x, "\\|"))),
      id = map_chr(header, 2) %>% str_remove(., ' gi'),
      prot_mystery_id = map_chr(header, 3),
      # seqs from genbank, refseq, emb, dbj
      prot_database = map_chr(header, 4),
      prot_accession = map_chr(header, 5),
      annotation = map_chr(header, 6),
      # which ICE
      parent_element = str_extract(annotation, "\\[.+\\]"),
      parent_element = str_remove_all(id, "\\[|\\]"),
      # protein annotation
      prot_name = str_remove(annotation, " \\[.+\\]")
    ) %>%
    select(-c(annotation, header)) %>%
    select(id, prot_name, everything())
}

prot_df <- MGEs %>% 
  mutate(protein_data = map(aa_fasta, parse_protein_fasta)) %>% 
  select(-contains('fasta')) %>% 
  unnest(protein_data)
prot_df


MGE_db <- full_join(dna_df, prot_df, by = c('type', 'id'))
glimpse(MGE_db)

write_delim(MGE_db, "./Data/ICE_db.tsv")

rm(dna_df, prot_df)

#