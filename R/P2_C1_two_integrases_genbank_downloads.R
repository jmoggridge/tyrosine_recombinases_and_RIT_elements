## P2_XX: Extract probable RIT elements using ID database

library(tidyverse)
library(progress)
library(glue)
library(rentrez)
source('./R/P2_parse_genbank_xml.R')
set_entrez_key('889fdb9786a14019a0a1257196a09ba4ba08')
Sys.getenv('ENTREZ_KEY')
dir.create('./data/CDD/genbank_2integrases/')

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
two_integrases <- 
  count_table |> 
  filter(n_prot_for_nuc ==2) |> 
  arrange(desc(n_prot_for_nuc)) |> 
  ungroup() |> 
  left_join(
    nuc_summary |> select(nuc_id, slen, organism, caption)
  )

# smallest nucleotides first; check which CD are present
two_integrases |>
  arrange(slen) |>  
  unnest(c(prot_id, cdd_title)) # |> View()

two_integrases |> unnest(c(prot_id, cdd_title)) |>  print(n=50)


write_rds(two_integrases, './data/CDD/ids_with_2_integrases.rds')

# how many nucleotides for a given taxa have > 3 RIT integrases?
nucs_by_taxid <- 
  two_integrases |> 
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
tot <- two_integrases |> 
  select(nuc_id) |> 
  distinct() |> 
  nrow()
tot

# skip these for now... problems parsing these genbank records
issues <- c()


## DOWNLOADS ------------------------------------------------------------------

## previously missing feature table
missing_cds <- read_rds('./data/CDD/rit_finder/missing_cds.rds') |> 
  select(nuc_id)

# fetch genbank records to see if CDS now present
pb <- make_pb(nrow(missing_cds))
missing_cds_fixed <- missing_cds |> 
  mutate(gbk = map(nuc_id, ~{
    print(.x)
    gb <- fetch_genbank_with_parts(.x)
    pb$tick()
    return(gb)
  }))
beep()
missing_cds_fixed
write_rds('data/CDD/RIT_gbk_noCDS_fixed.rds')

gbks <- read_rds('data/CDD/RIT_gbk_noCDS_fixed.rds')

# 2058 records with 2 CDD RIT (A,B,C) =============================

# 1,491 taxa not represented in rit elements already
to_retrieve1 <- 
  two_integrases |> 
  filter(!tax_id %in% c(read_rds('./results/rit_elements.rds') |> pull(tax_id))) |> 
  select(nuc_id, tax_id) |> 
  filter(!nuc_id %in% issues) |> 
  distinct()

# lower priority, taxa already found RITs in 
to_retrieve2 <- 
  two_integrases |> 
  filter(!nuc_id %in% to_retrieve1$nuc_id)

rm(cdd_titles, count_table, id_data, nucs_by_taxid, 
   nuc_summary, tot, two_integrases)


# processes a batch of requests and then saves the raw genbank data in the same row as id. No return value, just makes a file
fetch_gb_set <- function(df, id_var, start, stop, folder){
  pb <- make_pb(n = stop - start + 1)
  records <- df |> 
    slice(start:stop) |> 
    mutate(gbk = map({{nuc_id}}, ~{
      cat('\n')
      print(glue('id = .x'))
      gb <- fetch_genbank_with_parts(.x)
      pb$tick()
      return(gb)
    }))
  print(glimpse(records))
  write_rds(records, glue('{folder}/GB_w_parts_{start}_{stop}.rds'))
  beepr::beep()
}

# do fetching for sets of 100 ids.
purrr::walk(
  .x = seq(1, nrow(two_integrases), 100),
  .f = ~fetch_gb_set(
    df = two_integrases, 
    id_var = nuc_id, 
    start = .x, stop = .x + 99,
    folder = 'data/CDD/genbank_2integrases'
    )
)



