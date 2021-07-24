## P2_A4: Extract probable RIT elements using ID database

library(tidyverse)
library(progress)
library(glue)
library(rentrez)
source('./R/P2_parse_genbank_xml.R')
set_entrez_key('889fdb9786a14019a0a1257196a09ba4ba08')
Sys.getenv('ENTREZ_KEY')
dir.create('data/CDD/genbank')
#
id_data <- read_rds('./data/CDD/id_data_fixed.rds')
nuc_summary <- read_rds('./data/CDD/nuc_summary_fixed.rds')
cdd_titles <- read_rds('./data/CDD/cdd_summary.rds') |> 
  transmute(cdd_id, cdd_title = title)

# tally how many nucleotides/ taxa and how many proteins/ nucleotide
count_table <- 
  id_data |> 
  left_join(cdd_titles) |> 
  left_join(
    nuc_summary |> select(nuc_id, slen, organism, caption)
  ) |> 
  select(-cdd_id) |> 
  group_by(tax_id) |> 
  mutate(n_nuc_for_tax = length(nuc_id)) |> 
  group_by(tax_id, nuc_id) |> 
  summarize(
    n_prot_for_nuc = length(prot_id),
    prot_id = list(prot_id), 
    cdd_title = list(cdd_title)
    ) |> 
  ungroup()

count_table

# which nucleotides have > 3 RIT specific proteins?
three_integrases <- 
  count_table |> 
  filter(n_prot_for_nuc >=3) |> 
  arrange(desc(n_prot_for_nuc)) |> 
  ungroup() |> 
  left_join(
    nuc_summary |> select(nuc_id, slen, organism, caption)
    )

# smallest nucleotides first; check which CD are present
three_integrases |>
  arrange(slen) |>  
  unnest(c(prot_id, cdd_title)) # |> View()

three_integrases |> unnest(c(prot_id, cdd_title)) |>  print(n=50)


write_rds(three_integrases, 
          './data/CDD/ids_w_three_integrases.rds')

# how many nucleotides for a given taxa have > 3 RIT integrases?
nucs_by_taxid <- 
  three_integrases |> 
  group_by(tax_id) |> 
  summarise(
    n_nuc = length(nuc_id),
    nuc_id = list(nuc_id)
  ) |> 
  ungroup() |> 
  arrange(desc(n_nuc))
  
nucs_by_taxid |> 
  left_join(nuc_summary |> transmute(tax_id = taxid, organism) |> distinct())


# how many nucleotides with > 3 RIT specific proteins?
tot <- three_integrases |> 
  select(nuc_id) |> 
  distinct() |> 
  nrow()
tot

# skip these for now... problems parsing these genbank records
issues <-  c('1052702473', '483281088')
# three_integrases |> 
#   filter(nuc_id %in% issues)

pb <- progress_bar$new(total = 100)

to_retrieve <- 
  three_integrases |> 
  select(nuc_id, tax_id) |> 
  filter(!nuc_id %in% issues) |> 
  distinct()

# first 100
gb_records1 <- to_retrieve |> slice(1:100) |> 
  mutate(gbk = map(nuc_id, ~{
    pb$tick()
    fetch_genbank(.x)
    }))
beepr::beep()
write_rds(gb_records1, './data/CDD/genbank/RIT_gbk_1.rds')
rm(gb_records1)

# second 100
pb <- progress_bar$new(total = 100)
gb_records2 <- to_retrieve |> slice(101:200) |> 
  mutate(gbk = map(nuc_id, ~{
    pb$tick()
    fetch_genbank(.x)
  }))
beepr::beep()
write_rds(gb_records2, './data/CDD/genbank/RIT_gbk_2.rds')
rm(gb_records2)

# 3rd 100
pb <- progress_bar$new(total = 100)
gb_records3 <- to_retrieve |> slice(201:300) |> 
  mutate(gbk = map(nuc_id, ~{
    pb$tick()
    fetch_genbank(.x)
  }))
beepr::beep()
write_rds(gb_records3, './data/CDD/genbank/RIT_gbk_3.rds')
rm(gb_records3)

# 4th 100
pb <- progress_bar$new(total = 100)
gb_records4 <- to_retrieve |> slice(301:400) |> 
  mutate(gbk = map(nuc_id, ~{
    pb$tick()
    fetch_genbank(.x)
  }))
beepr::beep()
write_rds(gb_records4, './data/CDD/genbank/RIT_gbk_4.rds')
rm(gb_records4)

# 5th 100
pb <- progress_bar$new(total = 100)
gb_records5 <- to_retrieve |> slice(401:500) |> 
  mutate(gbk = map(nuc_id, ~{
    pb$tick()
    fetch_genbank(.x)
  }))
beepr::beep()
write_rds(gb_records5, './data/CDD/genbank/RIT_gbk_5.rds')
rm(gb_records5)

# last few records 100
pb <- progress_bar$new(total = 100)
gb_records6 <- to_retrieve |> slice(500:nrow(to_retrieve)) |> 
  mutate(gbk = map(nuc_id, ~{
    pb$tick()
    fetch_genbank(.x)
  }))
beepr::beep()
write_rds(gb_records6, './data/CDD/genbank/RIT_gbk_6.rds')
rm(gb_records6)

