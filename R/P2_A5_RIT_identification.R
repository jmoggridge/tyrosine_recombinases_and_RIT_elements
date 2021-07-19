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



create_errors <- c(
  # ids that seem to be causing parsing errors:
  # column `locus_tag` doesn't exist...
  "737980678", "47118329", "339284117", "1980667557", '47118328',
  # R just seems to hang when these are processed
  '56311475', '1008271296', '237624339', '338755570',
  '488453735', '227337254', '255767010', '184198282',
  '186474323', '1190193351', 
  
  '116222307'
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
  left_join(genbank_files_index)
  
# 531 records remaining for rit_finder
glimpse(three_ints_filter)

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

dir.create('./data/CDD/rit_finder/')
pb <- make_pb(200)
rit_list1 <- three_ints_filter |>
  dplyr::slice(1:200) |>
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
write_rds(rit_list1,
          './data/CDD/rit_finder/RIT_finder_rs_1_200.rds')



pb <- make_pb(50)
rit_list2 <- three_ints_filter |>
  dplyr::slice(201:400) |>
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

write_rds(rit_list2,
          './data/CDD/rit_finder/RIT_finder_rs_201_400.rds')
rm(rit_list2)

pb <- make_pb(n = nrow(three_ints_filter) - 400)
rit_list3 <- three_ints_filter |>
  dplyr::slice(401:nrow(three_ints_filter)) |>
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

write_rds(rit_list3, 
          './data/CDD/rit_finder/RIT_finder_rs_201_400.rds')
rm(rit_list3)

#
# TODO issue with many ids.... record is too large? need genbank records that actually contain CDS, some missing CDS
# seem to have cds when browsing but gives error with my code

