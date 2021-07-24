# P2_A5B: Rit_finder2 ----
# for "2 integrase" ids and previously missing data from P2_A4
# repeat process from P2_A5 and join to existing results....
# data/CDD/genbank/RIT_gbk_noCDS_fixed.rds
# and nucleotides with 2 integrases in CDD
# data/CDD/genbank_2integrases/GB_w_parts_....rds

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

### Run rit_finder ------
make_pb <- function(n){
  progress_bar$new(
    format = "Finding Rits: [:bar] :percent eta: :eta",
    total = n,
    clear = FALSE,
  )
}


## fixed records - CDS was missing before
fixed_cds <- read_rds('data/CDD/genbank/RIT_gbk_noCDS_fixed.rds') |> 
  unnest(gbk)
glimpse(fixed_cds)

# test:
# rit_finder2(fixed_cds$nuc_id[[1]], fixed_cds$gbk[[1]])

pb <- make_pb(nrow(fixed_cds))
rit_fixed_cds <- fixed_cds |>
  mutate(rit_output = map2(
    .x = nuc_id,
    .y = gbk,
    .f = ~{
      pb$tick()
      print(.x)
      rit_finder2(x = .x, gbk_xml = .y)
    }
  ))
rit_fixed_cds
rit_fixed_cds |> unnest()
beep()
write_rds(rit_fixed_cds, './data/CDD/rit_finder_rs/RIT_rs_for_fixed_cds.rds')



# these nucleotides probably contain RIT element bc have 2 CDD proteins.
two_ints <-
  read_rds('./data/CDD/ids_with_2_integrases.rds') |> 
  arrange(slen) |> 
  distinct()
glimpse(two_ints)

# create an index for nuc ids to their genbank records
genbank_files_index <- 
  tibble(file = Sys.glob('./data/CDD/genbank/GB_w_parts_*')) |> 
  mutate(nuc_id = map(file, ~ read_rds(.x) |> pull(nuc_id))) |> 
  unnest(nuc_id) |> 
  group_by(nuc_id) |> 
  filter(row_number() == 1) |> 
  ungroup()








