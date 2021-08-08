library(tidyverse)
library(rentrez)
library(glue)
library(progress)
library(tictoc)
library(beepr)

set_entrez_key('889fdb9786a14019a0a1257196a09ba4ba08')
Sys.getenv('ENTREZ_KEY')
# 


source('./R/P2_entrez_functions.R')
# MAIN: get ids ------------------------------------------------------

## CDD and prot ids first --------

# ## CDD search for list of RIT terms

terms <- tibble(term = map_chr(c('A', 'B', 'C'), ~glue('INT_Rit{.x}_C_like')))
cdd_searches <- terms |> cdd_search(term = term)
cdd_summary <- cdd_searches |> unnest_wider(cdd_summary)
write_rds(cdd_summary, './data/CDD/cdd_summary.rds')

## Link cdd records to proteins and nucleotides
cdd_prot_ids <-
  cdd_searches |>
  select(-cdd_summary) |>
  link_cdd_protein(id = cdd_id) |>
  unnest(prot_id)
beep()
write_rds(cdd_prot_ids, './data/CDD/cdd_prot_ids.rds')

cdd_prot_ids <- read_rds('./data/CDD/cdd_prot_ids.rds')


## Worked once, but usually fails...
cdd_prot_nuc_ids <-
  cdd_prot_ids |>
  link_protein_nuccore(id = prot_id) |>
  unnest(nuc_id)
beep()

write_rds(cdd_prot_nuc_ids, './data/CDD/cdd_prot_nuc_ids.rds')
rm(cdd_prot_ids)

cdd_prot_nuc_ids

id_data <-
  cdd_prot_nuc_ids |>
  link_nuccore_taxonomy(id = nuc_id)
beepr::beep()
write_rds(id_data, './data/CDD/cdd_prot_nuc_tax_ids.rds')
rm(cdd_prot_nuc_ids)


# Downloads ----
## check number of records before committing to downloads!!

id_data <- read_rds('./data/CDD/cdd_prot_nuc_tax_ids.rds')

## Taxonomy data ------

# retrieve taxonomy records for set of unique ids
tax_data <- 
  id_data |>
  select(tax_id) |> 
  filter(!is.na(tax_id)) |> 
  distinct() |> 
  fetch_taxonomy(id = tax_id)
print(tax_data)
beep()
write_rds(tax_data, './data/CDD/tax_data.rds')
rm(tax_data)  

## Prot data -------

# get protein data and summaries
prot_data <- id_data |>
  fetch_data(id = prot_id, db = 'protein', chunk_size = 50) |>
  unnest(cols = c(id, prot_name, prot_seq)) |>
  transmute(prot_id = id, prot_name, prot_seq)
beep()
write_rds(prot_data, './data/CDD/prot_data.rds')

prot_summary <- id_data |>
  fetch_summaries(id = prot_id, db = 'protein', chunk_size = 50)

prot_summary |> select(-token) |> unnest_wider(summary)
beep()
write_rds(prot_summary, './data/CDD/prot_summary.rds')
rm(prot_data, prot_summary, prot_ids)


## Nucleotide data  -----

nuc_ids <- id_data |>
  unnest(nuc_id) |> 
  select(nuc_id) |> 
  distinct()
rm(id_data)

dir.create('./data/CDD/nuc_data')

n <- nrow(nuc_ids)
nuc_sets <- 
  tibble(
    id_set = map(seq(1, n, 1000), ~nuc_ids[.x: min(.x + 999, n),])
    ) |> 
  mutate(set_number = row_number(),
         filepath = map2(
           .x = id_set, 
           .y = set_number,
           .f = ~fetch_nuc_datasets(data = .x, set_number = .y)
           )
         )


# uh oh some tibbles are NOT the right length, not enough seqs for ids.
# hard to unnest./join somehow?? Or redownload? to get missing.
# or maybe these seqs are removed so remove them from the id list...?

try(
  read_rds('./data/CDD/nuc_data/1.rds') |> 
    unnest(c(id, nuc_name, nuc_seq))
  )

# one chunk didn't work; missing a single sequence...
try1 <- 
  read_rds('./data/CDD/nuc_data/1.rds') |> 
  mutate(n_ids = map_dbl(id, length),
         n_seqs = map_dbl(nuc_seq, length)) |> 
  filter(n_ids == n_seqs) |> 
  select(-c(n_ids, n_seqs)) |> 
  unnest(c(id, nuc_name, nuc_seq))

print(try1)
rm(try1)

try1 |> mutate(kbp = nchar(nuc_seq)/1000) |> 
  ggplot(aes(kbp)) +
  geom_histogram() +
  scale_x_log10()

## TODO fix this issue in P2_A2 script because could be intensive.

## Nucleotide summaries ----

# unique ids
nuc_ids <- 
  id_data |> 
  select(nuc_id) |> 
  distinct()

# get summaries for first half using fetch_summaries
halfway <- nrow(nuc_ids) %/% 2

nuc_summary1 <- nuc_ids |>
  slice(1:halfway) |> 
  fetch_summaries(id = nuc_id, db = 'nuccore', chunk_size = 50) |> 
  select(-c(id, token)) |> 
  # force data to character vectors
  mutate(summary = map(summary, ~{
    .x |> mutate(across(everything(), as.character))
  })) |> 
  unnest(cols = c(summary))
beep()

# 2nd half
nuc_summary2 <- nuc_ids |>
  slice(halfway + 1: nrow(nuc_ids)) |> 
  fetch_summaries(id = nuc_id, db = 'nuccore', chunk_size = 50) |> 
  select(-c(id, token)) |> 
  mutate(summary = map(summary, ~{
    .x |> mutate(across(everything(), as.character))
  })) |> 
  unnest(cols = c(summary))
write_rds(nuc_summary1, './data/CDD/nuc_summary1.rds')

beep()

nuc_summary <- bind_rows(nuc_summary1, nuc_summary2)
write_rds(nuc_summary, './data/CDD/nuc_summary.rds')

rm(nuc_ids, nuc_summary, nuc_summary1, nuc_summary2, halfway)


# 
# ## for 16 individual sets....
# nuc_data1 <- nuc_ids |> slice(1:1000) |> 
#   fetch_data(id = nuc_id, db = 'nuccore', chunk_size = 100) |> 
#   parse_nuc_data()
# write_rds(nuc_data1, './data/CDD/nuc_data/1.rds', compress = 'gz')
# beep()
# rm(nuc_data1)
# 
# nuc_data2 <- nuc_ids |> slice(1001:2000) |> 
#   fetch_data(id = nuc_id, db = 'nuccore', chunk_size = 100) |> 
#   parse_nuc_data()
# write_rds(nuc_data2, './data/CDD/nuc_data/2.rds')
# beep()
# rm(nuc_data2)
# 
# nuc_data3 <- nuc_ids |> slice(2001:3000) |> 
#   fetch_data(id = nuc_id, db = 'nuccore', chunk_size = 100) |> 
#   parse_nuc_data()
# write_rds(nuc_data3, './data/CDD/nuc_data/3.rds')
# beep()
# rm(nuc_data3)
# 
# nuc_data4 <- nuc_ids |> slice(3001:4000) |> 
#   fetch_data(id = nuc_id, db = 'nuccore', chunk_size = 100)  |> 
#   parse_nuc_data()
# write_rds(nuc_data4, './data/CDD/nuc_data/4.rds')
# beep()
# rm(nuc_data4)
# 
# nuc_data5 <- nuc_ids |> slice(4001:5000) |> 
#   fetch_data(id = nuc_id, db = 'nuccore', chunk_size = 100) |> 
#   parse_nuc_data()
# write_rds(nuc_data5, './data/CDD/nuc_data/5.rds')
# beep()
# rm(nuc_data5)
# 
# nuc_data6 <- nuc_ids |> slice(5001:6000) |> 
#   fetch_data(id = nuc_id, db = 'nuccore', chunk_size = 100) |> 
#   parse_nuc_data()
# write_rds(nuc_data6, './data/CDD/nuc_data/6.rds')
# beep()
# rm(nuc_data6)
# 
# nuc_data7 <- nuc_ids |> slice(6001:7000) |> 
#   fetch_data(id = nuc_id, db = 'nuccore', chunk_size = 100) |> 
#   parse_nuc_data()
# write_rds(nuc_data7, './data/CDD/nuc_data/7.rds')
# rm(nuc_data7)
# beep()
# 
# nuc_data8 <- nuc_ids |> slice(7001:8000) |> 
#   fetch_data(id = nuc_id, db = 'nuccore', chunk_size = 100) |> 
#   parse_nuc_data()
# write_rds(nuc_data8, './data/CDD/nuc_data/8.rds')
# beep()
# rm(nuc_data8)
# 
# nuc_data9 <- nuc_ids |> slice(8001:9000) |> 
#   fetch_data(id = nuc_id, db = 'nuccore', chunk_size = 100) |> 
#   parse_nuc_data()
# write_rds(nuc_data9, './data/CDD/nuc_data/9.rds')
# rm(nuc_data9)
# beep()
# 
# nuc_data10 <- nuc_ids |> slice(9001:10000) |> 
#   fetch_data(id = nuc_id, db = 'nuccore', chunk_size = 100) |> 
#   parse_nuc_data()
# write_rds(nuc_data10, './data/CDD/nuc_data/10.rds')
# rm(nuc_data10)
# beep()
# 
# nuc_data11 <- nuc_ids |> slice(10001:11000) |> 
#   fetch_data(id = nuc_id, db = 'nuccore', chunk_size = 100) |> 
#   parse_nuc_data()
# write_rds(nuc_data11, './data/CDD/nuc_data/11.rds')
# rm(nuc_data11)
# beep()
# 
# nuc_data12 <- nuc_ids |> slice(11001:nrow(nuc_ids)) |> 
#   fetch_data(id = nuc_id, db = 'nuccore', chunk_size = 100) |> 
#   parse_nuc_data()
# beep()
# write_rds(nuc_data12, './data/CDD/nuc_data/12.rds')
# rm(nuc_data12)







