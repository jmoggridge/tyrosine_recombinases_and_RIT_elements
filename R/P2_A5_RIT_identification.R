## P2_A5 figure out how to identify rit element from genbank...

# uses functions from P2_rit_finder on the nucleotides with > 3 specific proteins (selected in P2_A3), for which genbank data was obtained in P2_A4.


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
dir.create('./data/CDD/rit_finder_rs/')

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
  tibble(file = Sys.glob('./data/CDD/genbank/RIT_gbk_[0-9].rds')) |> 
  mutate(nuc_id = map(file, ~ read_rds(.x) |> pull(nuc_id))) |> 
  unnest(nuc_id) |> 
  group_by(nuc_id) |> 
  filter(row_number() == 1) |> 
  ungroup()


have_genbank <- three_ints |> 
  filter(nuc_id %in% genbank_files_index$nuc_id)

# these are ids that genbank records aren't present for
no_genbank <- 
  anti_join(three_ints, have_genbank) |> pull(nuc_id) |> unique()
no_genbank


create_errors <- c(
  # ids that seem to be causing parsing errors
  # R just seems to hang when these are processed
  '56311475', '1008271296', '237624339', '338755570',
  '488453735', '227337254', '255767010', '184198282',
  '186474323', '1190193351', '116222307'
  )

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
  filter(!nuc_id %in% no_genbank) |> 
  left_join(genbank_files_index) |> 
  distinct()
  
# 531 records remaining for rit_finder
glimpse(three_ints_filter)
length(unique(three_ints_filter$nuc_id))

rm(no_genbank, have_genbank, three_ints, genbank_files_index)

## Main2 -----

### Run rit_finder ------
make_pb <- function(n){
  progress_bar$new(
    format = "Finding Rits: [:bar] :percent eta: :eta",
    total = n,
    clear = FALSE,
  )
}


pb <- make_pb(250)
rit_list1 <- three_ints_filter |>
  dplyr::slice(1:250) |>
  select(nuc_id, file) |>
  mutate(rit_output = map2(
    .x = nuc_id,
    .y = file,
    .f = ~{
      pb$tick()
      print(.x)
      rit_finder(x = .x, genbank_lib = .y)
      }
    ))
beep()
write_rds(rit_list1, './data/CDD/rit_finder_rs/RIT_finder_rs_1_250.rds')

pb <- make_pb(n = nrow(three_ints_filter) - 251 + 1)
rit_list2 <- three_ints_filter |>
  dplyr::slice(251:nrow(three_ints_filter)) |>
  select(nuc_id, file) |>
  mutate(rit_output = map2(
    .x = nuc_id,
    .y = file,
    .f = ~{
      pb$tick()
      print(.x)
      rit_finder(x = .x, genbank_lib = .y)
    }
  ))
beep()
write_rds(rit_list2, './data/CDD/rit_finder_rs/RIT_finder_rs_251_517.rds')


#
# TODO issue with many ids.... records too large? 
# need genbank records that actually contain CDS, some missing CDS
# seem to have cds when browsing but gives error with my code...FIXED in P@_A4B
# some have no column 'locus tag' which causes an error...FIXED
