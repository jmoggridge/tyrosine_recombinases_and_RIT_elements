## P2_A2_repair_nucleotide_data -----


library(tidyverse)
source('./R/P2_entrez_functions.R')

id_data <- read_rds('./data/CDD/cdd_prot_nuc_tax_ids.rds')
summaries <- read_rds('./data/CDD/nuc_summary.rds')
map(id_data, ~sum(is.na(.x)))

no_tax_id <- id_data |> 
  filter(is.na(tax_id))
  

nuc_ids <- 
  id_data |> 
  select(nuc_id) |> 
  filter(!is.na(nuc_id)) |> 
  distinct()

n <- nrow(nuc_ids)
nuc_sets <- 
  tibble(id_set = map(seq(1, n, 1000), ~nuc_ids[.x: min(.x + 999, n),])) |> 
  mutate(set_number = row_number())
nuc_sets

last <- read_rds('./data/CDD/nuc_data/12.rds') |> 
  select(id) |> 
  unnest(c(id))
last
