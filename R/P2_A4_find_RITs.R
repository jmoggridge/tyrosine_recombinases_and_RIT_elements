## P2_A4: Extract probable RIT elements using ID database

library(tidyverse)
library(progress)
source('./R/P2_parse_genbank_xml.R')

#
id_data <- read_rds('./data/CDD/id_data_fixed.rds')
nuc_summary <- read_rds('./data/CDD/nuc_summary_fixed.rds')

# tally how many nucleotides/ taxa and how many proteins/ nucleotide
count_table <- 
  id_data |> 
  select(-cdd_id) |> 
  group_by(tax_id) |> 
  mutate(n_nuc_for_tax = length(nuc_id)) |> 
  group_by(tax_id, nuc_id) |> 
  summarize(n_prot_for_nuc = length(prot_id),
         proteins = list(prot_id))
count_table
# which nucleotides have > 3 RIT specific proteins?
three_integrases <- count_table |> 
  filter(n_prot_for_nuc >=3) |> 
  arrange(desc(n_prot_for_nuc)) |> 
  ungroup()

three_integrases

# how many nucleotides with > 3 RIT specific proteins?
tot <- three_integrases |> 
  select(nuc_id) |> 
  distinct() |> 
  nrow()
tot

# skip these for now...
issues <-  c('1052702473', '483281088')
pb <- progress_bar$new(total = tot)

gb_records <- 
  three_integrases |> 
  filter(!nuc_id %in% issues) |> 
  select(nuc_id, tax_id) |> 
  distinct() |> 
  mutate(gbk = map(nuc_id, ~{
    pb$tick()
    fetch_genbank(.x)
    }))

write_rds(gb_records, './data/CDD/genbank_for_nuc_w_3_integrases')

# problems
# 1052702473
# 483281088
# 
