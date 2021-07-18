## P2_A5 figure out how to identify rit element from genbank...

## Setup ----

library(tidyverse)
library(progress)
library(glue)
library(Biostrings)
library(furrr)
library(tidymodels)
library(beepr)
library(tictoc)

# load functions
source('./R/P2_rit_finder.R')


## Main1 -----

### Load probable RIT elements id data ----

# genbank files for these nuc ids lack the CDS features.
no_cds <- c('1834417989', '1024364864', '1817592545')

# these nucleotides probably contain RIT element bc have 3 CDD proteins.
three_ints <- read_rds('./data/CDD/ids_w_three_integrases.rds') |> 
  filter(!nuc_id %in% no_cds) |> 
  arrange(slen)
glimpse(three_ints)

# create an index for nuc ids to their genbank records
genbank_files_index <- 
  tibble(file = Sys.glob('./data/CDD/RIT_gbk_[0-9].rds')) |> 
  mutate(nuc_id = map(file, ~ read_rds(.x) |> pull(nuc_id))) |> 
  unnest(nuc_id) |> 
  distinct()


## Main2 -----

### Run rit_finder ------

three_ints |> 
  filter(nuc_id %in% genbank)
# TODO issue with many ids.... record is too large? need genbank records that actually contain CDS, some missing CDS

# No cds issues with: 7, 11, 13, 14*, 16, 17, 19,22, 32*, 34, 35, 46... many others.
# seem to have cds when browsing but gives error with my code
creates_errors <- c(
  "737980678" ,# nuc_id [48],
  "47118329", # 219
  "339284117", # 
  "1980667557", # not in file index...
  ''
  )

# 
# pb <- progress_bar$new(
#   format = "Finding Rits: [:bar] :percent eta: :eta",
#   total = (nrow(three_ints) - length(creates_errors)),
#   clear = FALSE, 
#   )
# 
# rits_list_all_nuc <- three_ints |> 
#   filter(!nuc_id %in% creates_errors) |> 
#   pull(nuc_id) |> 
#   map(~{pb$tick()
#     rit_finder(.x)}) 
# beep()

length(rits_list)

rits_output <- 
  rits_list |> 
  purrr::reduce(bind_rows) |> 
  left_join(three_ints) |> 
  mutate(
    null_rits = map_lgl(rits, ~ifelse(is.null(nrow(.x)), T, F))
  ) |> 
  filter(null_rits == F) |> 
  mutate(n_rits = map_dbl(rits, nrow))

write_rds(rits_output, './data/CDD/Rits.rds')

# rits_output |>  View()
# rits_output |> unnest(rits) |> View()


# TODO figure out which ids are missing from the results
three_ints |> 
  select(nuc_id) |> 
  anti_join(rits_output |> select(nuc_id))





# pb <- progress_bar$new(total = 47)
# 
# rit_by_id_list1 <- three_ints |>
#   dplyr::slice(1:47) |>
#   pull(nuc_id) |>
#   map(~{pb$tick()
#       rit_finder(.x)})
# # beep()
# 
# 
# pb <- progress_bar$new(total = 52)
# 
# rit_by_id_list2 <- three_ints |>
#   dplyr::slice(49:100) |>
#   pull(nuc_id) |>
#   map(~{pb$tick()
#       rit_finder(.x)})
# # beep()
# # rit_by_id_list2
# 
# pb <- progress_bar$new(total = 100)
# rit_by_id_list3 <- three_ints |>
#   dplyr::slice(101:200) |>
#   pull(nuc_id) |>
#   map(~{pb$tick()
#     rit_finder(.x)})
# # beep()
# # rit_by_id_list3




# three_ints$nuc_id[219] "47118329"
# x Column `locus_tag` doesn't exist.
# "339284117" # 226
# Error in gzfile(file, "rb") : invalid 'description' argument


remaining <- three_ints |>
  dplyr::slice(201:300) |>
  filter(!nuc_id %in% creates_errors) |> 
  pull(nuc_id)

pb <- progress_bar$new(
  format = "Finding Rits: [:bar] :percent eta: :eta",
  total = length(remaining),
  clear = FALSE,
  )

rit_by_id_list4 <-
  remaining |> 
  map(~{
    pb$tick()
    print(.x)
    rit_finder(.x)
    })
beep()


pb <- progress_bar$new(total = 100)
rit_by_id_list5 <- three_ints |>
  dplyr::slice(301:400) |>
  pull(nuc_id) |>
  map(~{pb$tick()
    rit_finder(.x)})
# beep()
# rit_by_id_list5
# 

three_ints$nuc_id[410]
pb <- progress_bar$new(total = 100)
rit_by_id_list6 <- three_ints |>
  dplyr::slice(401:500) |>
  pull(nuc_id) |>
  map(~{pb$tick()
    rit_finder(.x)})
# 
# [======>---------------------------------------------------------------]  10%Error in gzfile(file, "rb") : invalid 'description' argument
# In addition: There were 11 warnings (use warnings() to see them)
beep()
# rit_by_id_list6












# # prot_id and accessions for this set of nt's.
# three_ints_prots <- three_ints |> unnest(prot_id) |> pull(prot_id)
# 
# # join with other accession numbers
# prot_accessions <- read_rds('./data/CDD/prot_summary_fixed.rds') |> 
#   filter(prot_id %in% three_ints_prots) |> 
#   transmute(prot_id, prot_accession = glue('{caption}.1'))
# 
# # join the protein accession labels that match the genbank CDS entries
# three_ints <- three_ints |> 
#   unnest(c(prot_id, cdd_title)) |> 
#   left_join(prot_accessions, by = c("prot_id")) |> 
#   group_by(nuc_id) |> 
#   nest(prot_data = c(prot_id, prot_accession, cdd_title)) |> 
#   mutate(nuc_accession = caption) |> 
#   ungroup()