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

# load functions and models
source('./R/P2_rit_finder.R')


## Main1 -----

### Load probable RIT elements id data ----

# these nucleotides probably contain RIT element bc have 3 CDD proteins.
# TODO go back to P2_A4 and find out why duplicated rows here
three_ints <- read_rds('./data/CDD/ids_w_three_integrases.rds') |> 
  arrange(slen) |> 
  distinct()
glimpse(three_ints)

# create an index for nuc ids to their genbank records
genbank_files_index <- 
  tibble(file = Sys.glob('./data/CDD/RIT_gbk_[0-9].rds')) |> 
  mutate(nuc_id = map(file, ~ read_rds(.x) |> pull(nuc_id))) |> 
  unnest(nuc_id) |> 
  distinct()

have_genbank <- three_ints |> 
  filter(nuc_id %in% genbank_files_index$nuc_id)

# these are ids that genbank records aren't present for
no_genbank <- 
  anti_join(three_ints, have_genbank) |> pull(nuc_id) |> unique()
no_genbank


# ids that are causing errors
create_errors <- c("737980678", "47118329", "339284117")

## ERRORS
##  "47118329"
##  x Column `locus_tag` doesn't exist.
##  "339284117", "737980678"
##  Error in gzfile(file, "rb") : invalid 'description' argument

# genbank files for these nuc ids lack the CDS features.
# no_CDS <- c('1834417989', '1024364864', '1817592545')
# No cds issues with many records: 7, 11, 13, 14*, 16, 17, 19,22, 32*, 34, 35, 46... many others.


# remove any ids that have no genbank or cause other errors
three_ints_filter <- three_ints |> 
  filter(!nuc_id %in% create_errors) |> 
  left_join(genbank_files_index)
  
# 531 records remaining for rit_finder
glimpse(three_ints_filter)


## Main2 -----

### Run rit_finder ------
make_pb <- function(n){
  progress_bar$new(
    format = "Finding Rits: [:bar] :percent eta: :eta",
    total = n,
    clear = FALSE,
  )
}

pb <- make_pb(100)

rit_list1 <- three_ints_filter |> 
  dplyr::slice(1:100) |> 
  select(nuc_id, file) |> 
  mutate(rit_output = map2(
    .x = nuc_id, 
    .y = file, 
    .f = ~{
      pb$tick()
      print(.x)
      rit_finder(x = .x, genbank_library = .y)
      }
    ))
beep()

pb <- make_pb(100)

rit_list2 <- three_ints_filter |> 
  dplyr::slice(101:200) |> 
  select(nuc_id, file) |> 
  mutate(rit_output = map2(
    .x = nuc_id, 
    .y = file, 
    .f = ~{
      pb$tick()
      print(.x)
      rit_finder(x = .x, genbank_library = .y)
    }
  ))


pb <- make_pb(100)
rit_list3 <- three_ints_filter |> 
  dplyr::slice(201:300) |> 
  select(nuc_id, file) |> 
  mutate(rit_output = map2(
    .x = nuc_id, 
    .y = file, 
    .f = ~{
      pb$tick()
      print(.x)
      rit_finder(x = .x, genbank_library = .y)
    }
  ))

pb <- make_pb(100)

rit_list4 <- three_ints_filter |> 
  dplyr::slice(301:400) |> 
  select(nuc_id, file) |> 
  mutate(rit_output = map2(
    .x = nuc_id, 
    .y = file, 
    .f = ~{
      pb$tick()
      print(.x)
      rit_finder(x = .x, genbank_library = .y)
    }
  ))

pb <- make_pb(100)

rit_list5 <- three_ints_filter |> 
  dplyr::slice(401:500) |> 
  select(nuc_id, file) |> 
  mutate(rit_output = map2(
    .x = nuc_id, 
    .y = file, 
    .f = ~{
      pb$tick()
      print(.x)
      rit_finder(x = .x, genbank_library = .y)
    }
  ))

# TODO issue with many ids.... record is too large? need genbank records that actually contain CDS, some missing CDS
# seem to have cds when browsing but gives error with my code





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




# 
# remaining <- three_ints |>
#   dplyr::slice(201:300) |>
#   filter(!nuc_id %in% creates_errors) |> 
#   pull(nuc_id)
# 
# pb <- progress_bar$new(
#   format = "Finding Rits: [:bar] :percent eta: :eta",
#   total = length(remaining),
#   clear = FALSE,
#   )
# 
# rit_by_id_list4 <-
#   remaining |> 
#   map(~{
#     pb$tick()
#     print(.x)
#     rit_finder(.x)
#     })
# beep()
# 
# 
# pb <- progress_bar$new(total = 100)
# rit_by_id_list5 <- three_ints |>
#   dplyr::slice(301:400) |>
#   pull(nuc_id) |>
#   map(~{pb$tick()
#     rit_finder(.x)})
# # beep()
# # rit_by_id_list5
# # 
# 
# three_ints$nuc_id[410]
# pb <- progress_bar$new(total = 100)
# rit_by_id_list6 <- three_ints |>
#   dplyr::slice(401:500) |>
#   pull(nuc_id) |>
#   map(~{pb$tick()
#     rit_finder(.x)})
# # 
beep()
# rit_by_id_list6


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
# 
# length(rits_list)
# 
# rits_output <- 
#   rits_list |> 
#   purrr::reduce(bind_rows) |> 
#   left_join(three_ints) |> 
#   mutate(
#     null_rits = map_lgl(rits, ~ifelse(is.null(nrow(.x)), T, F))
#   ) |> 
#   filter(null_rits == F) |> 
#   mutate(n_rits = map_dbl(rits, nrow))
# 
# write_rds(rits_output, './data/CDD/Rits.rds')
# 
# # TODO figure out which ids are missing from the results
# three_ints |> 
#   select(nuc_id) |> 
#   anti_join(rits_output |> select(nuc_id))

