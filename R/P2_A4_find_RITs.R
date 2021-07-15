## P2_A4: Extract probable RIT elements using ID database

library(tidyverse)
library(progress)
library(rentrez)
source('./R/P2_parse_genbank_xml.R')
set_entrez_key('889fdb9786a14019a0a1257196a09ba4ba08')
Sys.getenv('ENTREZ_KEY')

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

# skip these for now... problems parsing these genbank records
issues <-  c('1052702473', '483281088')
pb <- progress_bar$new(total = 100)

to_retrieve <- 
  three_integrases |> 
  select(nuc_id, tax_id) |> 
  filter(!nuc_id %in% issues) |> 
  distinct()

# first 100
gb_records1 <- 
  to_retrieve |> 
  slice(1:100) |> 
  mutate(gbk = map(nuc_id, ~{
    pb$tick()
    fetch_genbank(.x)
    }))
beepr::beep()
write_rds(gb_records1, './data/CDD/RIT_gbk_1.rds')

# second 100
gb_records2 <- 
  to_retrieve |> 
  slice(101:200) |> 
  mutate(gbk = map(nuc_id, ~{
    pb$tick()
    fetch_genbank(.x)
  }))
beepr::beep()
write_rds(gb_records2, './data/CDD/RIT_gbk_2.rds')

# 3rd 100
gb_records3 <- 
  to_retrieve |> 
  slice(201:300) |> 
  mutate(gbk = map(nuc_id, ~{
    pb$tick()
    fetch_genbank(.x)
  }))
beepr::beep()
write_rds(gb_records3, './data/CDD/RIT_gbk_3.rds')

# 4th 100
gb_records4 <- 
  to_retrieve |> 
  slice(301:400) |> 
  mutate(gbk = map(nuc_id, ~{
    pb$tick()
    fetch_genbank(.x)
  }))
beepr::beep()
write_rds(gb_records4, './data/CDD/RIT_gbk_4.rds')

# 5th 100
gb_records5 <- 
  to_retrieve |> 
  slice(401:500) |> 
  mutate(gbk = map(nuc_id, ~{
    pb$tick()
    fetch_genbank(.x)
  }))
beepr::beep()
write_rds(gb_records5, './data/CDD/RIT_gbk_5.rds')
