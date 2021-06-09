

library(tidyverse)
library(Biostrings)
library(rhmmer)

# All the full protein sequences for 20 subfamilies from SMART
integrases <- 
  tibble(path = Sys.glob('./data/SMART/full_protein_fasta/*.fasta')) |> 
  mutate(
    domain_name = str_remove_all(path, './data/SMART/full_protein_fasta/|\\_proteins.fasta'),
    fasta = map(path, readAAStringSet),
    n_protein = map_int(fasta, ~length(.x)),
    prot_name = map(fasta, ~names(.x)),
    prot_seq = map(fasta, paste)
  ) |> 
  # add accession numbers
  left_join(read_csv('./data/SMART/SMART_YR_ids.csv') |> 
              janitor::clean_names(), by = 'domain_name',
  ) 

write_rds(integrases, './data/SMART_integrases.rds')

integrases |>  
  select(domain_name, id, prot_name) |> 
  unnest(cols = prot_name) |> 
  mutate(prot)

integrases



read_hmmersearch_tbl <- function(path){
  read_table(
    file = path,
    na = '-',
    comment = '#',
    col_names =  c(
      "prot_name", "prot_accession", 
      "hmm_name", "hmm_accesion",
      "seq_eval", "seq_score", "seq_bias",
      "best_dom_eval", "best_dom_score", "best_dom_bias",
      'domain_n_exp', 'domain_n_reg', 'domain_n_clu', 'domain_n_ov',
      'domain_n_env', 'domain_n_dom', 'domain_n_rep', 'domain_n_inc',
      "description")
  ) |> 
    select(prot_name, hmm_name, best_dom_score) |> 
    distinct() |> 
    pivot_wider(names_from = hmm_name, values_from = best_dom_score) |> 
    unnest(cols = c())
}


# one table
subfam <- 'Arch1'

eg_tab <- 
  read_hmmersearch_tbl("./data/SMART/hmmsearch_res/Arch1.search.tbl")
  
eg_tab

left_join()






# all hmmsearch output
hmmsearch_rs <- 
  tibble(path = Sys.glob('./data/SMART/hmmsearch_res/*.tbl')) |> 
  mutate(domain_name = str_extract(path, '[^\\s|!/]+.search.tbl'),
         domain_name = str_remove(domain_name, '.search.tbl')) |> 
  select(domain_name, path) |> 
  arrange(domain_name) |> 
  mutate(result = map(path, ~read_hmmertbl(.x)))
  # unnest(cols = c(result))

glimpse(hmmsearch_rs)


search_tbl <-
  rhmmer::read_tblout("./data/SMART/hmmsearch_res/Arch1.search.tbl") |> 
  select(-contains('accession'), -contains('domain_number_')) |>
  mutate(protein = paste(domain_name, description))

names(search_tbl)
glimpse(search_tbl)


