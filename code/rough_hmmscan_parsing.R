library(tidyverse)
library(Biostrings)
library(rhmmer)

# domtblout <- rhmmer::read_domtblout('./data/SMART/hmmscan_results/smart.scan.domtbl')
# glimpse(domtblout)

tblout <- rhmmer::read_tblout('./data/SMART/hmmscan_results/smart.scan.tbl')
glimpse(tblout)

# drop empty columns
tblout <- tblout |> 
    select(-contains('accession'))
glimpse(tblout)

# query names are clipped at the first space
# tblout |> 
#   pull(query_name) |> 
#   head()

# # eda
# tblout |> 
#   ggplot(aes(x=fct_rev(domain_name))) +
#   geom_bar() + coord_flip()
# tblout |> 
#   ggplot(aes(y = fct_rev(domain_name), x = best_domain_score)) +
#   geom_violin()
#   
#   

# All the full protein sequences for 20 subfamilies from SMART
smart_db <- 
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
  )  |>  
  select(domain_name, id, contains('n_'), everything()) 

smart_db

# make df with protein name key to match with hmmscan tblout
smart_prots <- smart_db |> 
  select(domain_name, n_protein, prot_name) |> 
  unnest(cols = c(prot_name)) |> 
  # take fasta header up until first space to match hmmscan
  mutate(query_name = trimws(str_extract(prot_name, '.*? '))) 

glimpse(smart_prots)


## join the hmmscan results to the proteins from smart
glimpse(tblout)

hmmscan_rs <- 
  left_join(tblout, smart_prots, by = 'query_name') |> 
  sample_n(1000)

pvt <- hmmscan_rs |> 
  arrange(domain_name) |> 
  select(domain_name, query_name, best_domain_score) |> 
  pivot_wider(names_from = domain_name, values_from = best_domain_score) |> 
  group_by(query_name) |> 
  summarize(across(Arch1:Xer, ~max(.x)))

glimpse(pvt)




