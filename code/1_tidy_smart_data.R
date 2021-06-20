library(Biostrings)
library(DECIPHER)
library(tidyverse)
library(janitor)

# make dataframe of fastas for each subfamily of domains
smart_domains <- 
  tibble(dom_path = Sys.glob('./data/SMART/domain_fasta/*.fasta')) |>
  mutate(
    subfamily = str_remove_all(dom_path, './data/SMART/domain_fasta/|\\.fasta'),
    dom_fasta = map(dom_path, readAAStringSet),
    dom_n = map_int(dom_fasta, ~length(.x)),
    dom_name = map(dom_fasta, ~names(.x)),
    dom_seq = map(dom_fasta, paste)
  ) |> 
  # add accession numbers
  left_join(read_csv('./data/SMART/smart_subfamily_ids.csv',
                     col_types = cols()) |> 
              janitor::clean_names() |> 
              select(-n_domains),
            by = 'subfamily') |> 
  select(subfamily, id, everything())

# All the full protein sequences for 20 subfamilies from SMART
smart_integrase_proteins <- 
  tibble(
    prot_path = Sys.glob('./data/SMART/full_protein_fasta/*.fasta')) |> 
  mutate(
    subfamily = str_remove_all(
      prot_path, './data/SMART/full_protein_fasta/|\\_proteins.fasta'),
    prot_fasta = map(prot_path, readAAStringSet),
    prot_n = map_int(prot_fasta, ~length(.x)),
    prot_name = map(prot_fasta, ~names(.x)),
    prot_seq = map(prot_fasta, paste)
  ) |> 
  # add accession numbers
  left_join(read_csv('./data/SMART/smart_subfamily_ids.csv', 
                     col_types = cols()) |> 
              janitor::clean_names() |> 
              select(-n_domains),
            by = 'subfamily') |> 
  select(subfamily, id, everything())

# join the full proteins and domains data together by subfamily
smart_df <- 
  left_join(smart_domains, 
            smart_integrase_proteins, by = c('subfamily', 'id')) |> 
  mutate(subfamily = as_factor(subfamily)) |> 
  select(-contains('path'))

glimpse(smart_df)

write_rds(smart_df, './data/smart_df.rds', compress = 'gz')




