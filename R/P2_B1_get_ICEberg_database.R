library(Biostrings)
library(tidyverse)
 
## The ICEberg2.0 database has a set of proteins from integrative conjugative genetic elements (ICEs) and other MGEs. https://db-mml.sjtu.edu.cn/ICEberg2/download.html. Here, I obtained  data for all experimental-validated elements with intact sequences in DNA and protein fastas for ICEs, IMEs, AICEs, CIMEs, T4SS-type ICEs. There are a total of 16,752 proteins (predicted at least) in the full dataset, from 
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


# Extract ICE seqs and header labels from fasta (ea. row is an element)
parse_dna_fasta <- function(dna_fasta){
  df <- tibble(
    dna_name = names(dna_fasta),
    dna_seq = paste(dna_fasta)
  ) |> 
    separate(
      dna_name, 
      into = c(NA, 'id', 'name', 'database', 'accession',
               'description'), 
      sep = '\\|'
    ) |> 
    mutate(
      start_pos = str_extract(description, "^[0-9]+"),
      end_pos = str_extract(description, "\\.\\.[0-9|]+"),
      end_pos = str_remove_all(end_pos, '\\.'),
    ) |>
    select(name, id, database, accession, start_pos,
           end_pos, description, dna_seq)
}
# Extract protein seqs and header labels from fasta (multiple/element)
parse_protein_fasta <- function(aa_fasta){
  df <- 
    tibble(
      header = names(aa_fasta),
      prot_seq = paste(aa_fasta)
    ) |>
    separate(
      header, sep = '\\|', remove = T, 
      into = c(NA, 'id', 'prot_mystery_id', 'prot_database',
               'prot_accession', 'annotation'), 
    ) |> 
    #   # extract annotation info from sequence headers
    mutate(
      annotation = str_trim(annotation),
      # which ICE
      parent_element = str_extract(annotation, "\\[.+\\]"),
      parent_id = str_remove_all(id, "\\sgi"),
      # protein annotation
      prot_name = str_remove(annotation, " \\[.+\\]")
    ) |>
    select(-c(annotation, id)) |>
    select(parent_id, prot_name, everything())
}


## MAIN =============================================


# Organize data for all types of elements in single dataframe
MGEs <- 
  tibble(type = c('AICE', 'CIME', 'ICE', 'IME', 'T4SS-type_ICE')) |> 
  mutate(
    dna_fasta = map(type, ~get_ICEberg(.x, "seq")),
    aa_fasta = map(type, ~get_ICEberg(.x, "aa"))
    )

# check that all sequences are obtained
MGEs |> 
  mutate(n_dna = map_int(dna_fasta, length),
         n_protein = map_int(aa_fasta, length))
# # A tibble: 5 x 5
# type          dna_fasta  aa_fasta   n_dna n_protein
# <chr>         <list>     <list>     <int>     <int>
# 1 AICE          <DNAStrnS> <AAStrngS>     7       113
# 2 CIME          <DNAStrnS> <AAStrngS>     3        40
# 3 ICE           <DNAStrnS> <AAStrngS>   113      7493
# 4 IME           <DNAStrnS> <AAStrngS>    23       624
# 5 T4SS-type_ICE <DNAStrnS> <AAStrngS>   100      6555
# 
# Check total number of elements and proteins
MGEs |>
  mutate(n_dna = map_int(dna_fasta, length),
         n_protein = map_int(aa_fasta, length)) |> 
  summarise(elements = sum(n_dna),
            proteins = sum(n_protein))
# A tibble: 1 x 2
# elements proteins
# <int>    <int>
#   1      246    14825

# save raw nested dataset to file
write_rds(MGEs, "./data/iceberg/ICEberg_experimental_intact_database_raw.rds")


dna_df <- MGEs |> 
  mutate(dna_data = map(dna_fasta, parse_dna_fasta)) |> 
  select(-contains('fasta')) |> 
  unnest(dna_data)


dna_df |>
  ggplot(aes(y = type)) +
  geom_bar()

prot_df <- MGEs |> 
  mutate(protein_data = map(aa_fasta, parse_protein_fasta)) |> 
  select(-contains('fasta')) |> 
  unnest(protein_data) |> 
  mutate(id = parent_id)
prot_df

# combine the full elements (dna fasta data) with their proteins
MGE_db <- full_join(dna_df, prot_df, by = c('type', 'id'))

glimpse(MGE_db)


write_tsv(MGE_db, "./data/iceberg/ICE_db.tsv")

rm(dna_df, prot_df)

# # Fasta files
# dna_fasta_db <- purrr::reduce(MGEs$dna_fasta, ~c(.x,.y))
# length(dna_fasta_db)
# aa_fasta_db <- purrr::reduce(MGEs$aa_fasta, ~c(.x,.y))
# length(aa_fasta_db)
# 
# writeXStringSet(dna_fasta_db, filepath = './data/iceberg/iceberg_dna.fa')
# writeXStringSet(aa_fasta_db, filepath = './data/iceberg/iceberg_aa.fa')
# 






