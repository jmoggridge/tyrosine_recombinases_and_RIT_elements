# P2_A5B: Rit_finder2 ----

# for "2 integrase" ids and previously missing data from P2_A4
# repeat process from P2_A5 and join to existing results....
# data/CDD/genbank/RIT_gbk_noCDS_fixed.rds
# and nucleotides with 2 integrases in CDD
# data/CDD/genbank_2integrases/GB_w_parts_....rds

library(Biostrings)
library(tidyverse)
library(progress)
library(glue)
library(furrr)
library(tidymodels)
library(beepr)
library(tictoc)

# load parse_genbank function
source('./R/P2_parse_genbank_xml.R')
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

# genbank files too large to read xml. Need to fix parser to handle xml > 10Mb
cause_error <- c(
  "391353422", "1207852062", '58012118', '93352797', '93352797', '94308945',
  '120326872', '120537037'
  )

## fixed records - CDS was missing before
fixed_cds <- read_rds('data/CDD/genbank/RIT_gbk_noCDS_fixed.rds') |> 
  unnest(gbk) |> 
  mutate(size= map_dbl(gbk, ~object.size(.x)))
glimpse(fixed_cds)

# small genbank files
small_fixed_cds <- fixed_cds |> 
  filter(size < 1e7)
small_fixed_cds |> 
  print(n=100)

# medium sized ones
medium_fixed_cds <- fixed_cds |> 
  anti_join(small_fixed_cds,  by = c("nuc_id", "gbk", "size")) |> 
  filter(size < 1.7e7)

# large ones

## these ones aren't parsed by xml2 properly because they are too big! > 10Mb.
# going to try setting "HUGE" option on parse xml..
large_fixed_cds <- fixed_cds |> 
  anti_join(small_fixed_cds,  by = c("nuc_id", "gbk", "size")) |> 
  filter(size > 1.7e7)

rm(fixed_cds)


# small
pb <- make_pb(nrow(small_fixed_cds))
rit_fixed_cds_small <- small_fixed_cds |>
  mutate(rit_output = map2(
    .x = nuc_id,
    .y = gbk,
    .f = ~{
      pb$tick()
      print(.x)
      rit_finder2(x = .x, gbk_xml = .y)
    }
  ))
rit_fixed_cds_small
rit_fixed_cds_small |> unnest() |> unnest()
beep()
write_rds(rit_fixed_cds_small, './data/CDD/rit_finder_rs/RIT_rs_for_fixed_cds_small.rds')
rm(rit_fixed_cds_small, small_fixed_cds)

# medium
pb <- make_pb(nrow(medium_fixed_cds))
rit_fixed_cds_med <- medium_fixed_cds |>
  mutate(rit_output = map2(
    .x = nuc_id,
    .y = gbk,
    .f = ~{
      pb$tick()
      print(.x)
      rit_finder2(x = .x, gbk_xml = .y)
    }
  ))
beep()
rit_fixed_cds_med
rit_fixed_cds_med |> unnest() |> unnest()
beep()
write_rds(rit_fixed_cds_med,
          './data/CDD/rit_finder_rs/RIT_rs_for_fixed_cds_med.rds')
rm(rit_fixed_cds_med, medium_fixed_cds)


#### large sized ones - not working genbank has different format when parsing huge file - newline characters are creating unnecessary columns.
### # need to fix parser to work with read_xml(..., option = 'HUGE')
# 
# gbk <- large_fixed_cds$gbk[[1]]
# # x <- parse_genbank(gbk = gbk)
# gbk_xml <- read_xml(gbk, options = 'HUGE') |> as_list()
# gbk <- gbk_xml
# 
# pb <- make_pb(5)
# rit_fixed_cds_lrg <- large_fixed_cds |>
#   dplyr::slice(1:5) |>
#   mutate(rit_output = map2(
#     .x = nuc_id,
#     .y = gbk,
#     .f = ~{
#       pb$tick()
#       print(.x)
#       rit_finder2(x = .x, gbk_xml = .y)
#     }
#   ))
# beep()
# rit_fixed_cds_lrg |> unnest()
# rit_fixed_cds_lrg |> unnest() |> unnest()
# beep()
# write_rds(rit_fixed_cds_med,
#           './data/CDD/rit_finder_rs/RIT_rs_for_fixed_cds_lrg.rds')
# 


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

rm(large_fixed_cds, small_fixed_cds, medium_fixed_cds)

## Part 2 -------

# using rit_finder2 to extract trios from all the genbank records from the 'two integrases' set of nucleotide ids
# still have issues with some genbank >10 Mb. Some others contain no Rit elements.

##1 ----
gb_with_parts1 <- read_rds('./data/CDD/genbank/GB_w_parts_1_100.rds') |>
  unnest(gbk)

pb <- make_pb(nrow(gb_with_parts1))
rits_2int_1 <- gb_with_parts1 |>
  mutate(rit_output = map2(
    .x = nuc_id,
    .y = gbk,
    .f = ~{
      pb$tick()
      print(.x)
      rit_finder2(x = .x, gbk_xml = .y)
    }
  ))
beep()
rits_2int_1
rits_2int_1 |> unnest() |> unnest()
beep()
write_rds(rits_2int_1,
          './data/CDD/rit_finder_rs/rits_2int_1.rds')


## 2 ----

gb_with_parts2 <- read_rds('./data/CDD/genbank/GB_w_parts_101_200.rds') |>
  filter(!nuc_id %in% cause_error) |>
  unnest(gbk)

pb <- make_pb(nrow(gb_with_parts2))
rits_2int_2 <- gb_with_parts2 |>
  mutate(rit_output = map2(
    .x = nuc_id, .y = gbk,
    .f = ~{
      pb$tick()
      print(.x)
      rit_finder2(x = .x, gbk_xml = .y)
    }
  ))
beep()
rits_2int_2
rits_2int_2 |> unnest() |> unnest()
beep()
write_rds(rits_2int_2,
          './data/CDD/rit_finder_rs/rits_2int_2.rds')
rm(gb_with_parts2,
   rits_2int_2)

## 3 ----
gb_with_parts3 <- read_rds('./data/CDD/genbank/GB_w_parts_201_300.rds') |>
  filter(!nuc_id %in% cause_error) |>
  unnest(gbk) |>
  mutate(size= map_dbl(gbk, ~object.size(.x))) |>
  filter(size < 1.7e7)


pb <- make_pb(nrow(gb_with_parts3))
rits_2int_3 <- gb_with_parts3 |>
  mutate(rit_output = map2(
    .x = nuc_id, .y = gbk,
    .f = ~{
      pb$tick()
      print(.x)
      rit_finder2(x = .x, gbk_xml = .y)
    }
  ))
beep()
rits_2int_3
rits_2int_3 |> unnest() |> unnest()
beep()
write_rds(rits_2int_3,
          './data/CDD/rit_finder_rs/rits_2int_3.rds')
rm(gb_with_parts3,
   rits_2int_3)


## 4 -----
gb_with_parts4 <- read_rds('./data/CDD/genbank/GB_w_parts_301_400.rds') |>
  unnest(gbk) |>
  filter(!nuc_id %in% cause_error) |>
  mutate(size= map_dbl(gbk, ~object.size(.x))) |>
  filter(size < 1.7e7)

pb <- make_pb(nrow(gb_with_parts4))
rits_2int_4 <- gb_with_parts4 |>
  mutate(rit_output = map2(
    .x = nuc_id, .y = gbk,
    .f = ~{
      pb$tick()
      print(.x)
      rit_finder2(x = .x, gbk_xml = .y)
    }
  ))
beep()
rits_2int_4
rits_2int_4 |> unnest() |> unnest()
beep()
write_rds(rits_2int_4, './data/CDD/rit_finder_rs/rits_2int_4.rds')
rm(rits_2int_4)


## 5 -----
gb_with_parts5 <- read_rds('./data/CDD/genbank/GB_w_parts_401_500.rds') |>
  unnest(gbk) |>
  filter(!nuc_id %in% cause_error) |>
  mutate(size= map_dbl(gbk, ~object.size(.x))) |>
  filter(size < 1.7e7)


pb <- make_pb(nrow(gb_with_parts5))
rits_2int_5 <- gb_with_parts5 |>
  mutate(rit_output = map2(
    .x = nuc_id, .y = gbk,
    .f = ~{
      pb$tick()
      print(.x)
      rit_finder2(x = .x, gbk_xml = .y)
    }
  ))
beep()
rits_2int_5
rits_2int_5 |> unnest() |> unnest()
beep()
write_rds(rits_2int_5, './data/CDD/rit_finder_rs/rits_2int_5.rds')
rm(rits_2int_5, gb_with_parts5)



# 6 ------
gb_with_parts6 <- read_rds('./data/CDD/genbank/GB_w_parts_501_600.rds') |>
  unnest(gbk) |>
  filter(!nuc_id %in% cause_error) |>
  mutate(size= map_dbl(gbk, ~object.size(.x))) |>
  filter(size < 1.7e7)


pb <- make_pb(nrow(gb_with_parts6))
rits_2int_6 <- gb_with_parts6 |>
  mutate(rit_output = map2(
    .x = nuc_id, .y = gbk,
    .f = ~{
      pb$tick()
      print(.x)
      rit_finder2(x = .x, gbk_xml = .y)
    }
  ))
beep()
rits_2int_6
rits_2int_6 |> unnest() |> unnest()
beep()
write_rds(rits_2int_6, './data/CDD/rit_finder_rs/rits_2int_6.rds')
rm(rits_2int_6, gb_with_parts6)

## 7 ------
gb_with_parts7 <- read_rds('./data/CDD/genbank/GB_w_parts_601_700.rds') |>
  unnest(gbk) |>
  filter(!nuc_id %in% cause_error) |>
  mutate(size= map_dbl(gbk, ~object.size(.x))) |>
  filter(size < 1.7e7)


pb <- make_pb(nrow(gb_with_parts7))
rits_2int_7 <- gb_with_parts7 |>
  mutate(rit_output = map2(
    .x = nuc_id, .y = gbk,
    .f = ~{
      pb$tick()
      print(.x)
      rit_finder2(x = .x, gbk_xml = .y)
    }
  ))
beep()
rits_2int_7
rits_2int_7 |> unnest() |> unnest()
beep()
write_rds(rits_2int_7, './data/CDD/rit_finder_rs/rits_2int_7.rds')
rm(rits_2int_7, gb_with_parts7)


## 8 ----
gb_with_parts8 <-
  read_rds('./data/CDD/genbank/GB_w_parts_701_800.rds') |>
  unnest(gbk) |>
  filter(!nuc_id %in% cause_error) |>
  mutate(size= map_dbl(gbk, ~object.size(.x))) |>
  filter(size < 1.7e7)

pb <- make_pb(nrow(gb_with_parts8))
rits_2int_8 <- gb_with_parts8 |>
  mutate(rit_output = map2(
    .x = nuc_id, .y = gbk,
    .f = ~{
      pb$tick()
      print(.x)
      rit_finder2(x = .x, gbk_xml = .y)
    }
  ))
beep()
rits_2int_8
rits_2int_8 |> unnest() |> unnest()
beep()
write_rds(rits_2int_8, './data/CDD/rit_finder_rs/rits_2int_8.rds')
rm(rits_2int_8, gb_with_parts8)


## 9 ----

gb_with_parts9 <- read_rds('./data/CDD/genbank/GB_w_parts_801_900.rds') |>
  unnest(gbk) |> 
  filter(!nuc_id %in% cause_error) |>
  mutate(size= map_dbl(gbk, ~object.size(.x))) |> 
  filter(size < 1.7e7)


pb <- make_pb(nrow(gb_with_parts9))
rits_2int_9 <- gb_with_parts9 |>
  mutate(rit_output = map2(
    .x = nuc_id, .y = gbk,
    .f = ~{
      pb$tick()
      print(.x)
      rit_finder2(x = .x, gbk_xml = .y)
    }
  ))
beep()
rits_2int_9
rits_2int_9 |> unnest() |> unnest()
beep()
write_rds(rits_2int_9, './data/CDD/rit_finder_rs/rits_2int_9.rds')
rm(rits_2int_9, gb_with_parts9)


## 10 ----

gb_with_parts10 <- read_rds('./data/CDD/genbank/GB_w_parts_901_1000.rds') |>
  unnest(gbk) |> 
  filter(!nuc_id %in% cause_error) |>
  mutate(size= map_dbl(gbk, ~object.size(.x))) |> 
  filter(size < 1.7e7)

pb <- make_pb(nrow(gb_with_parts10))
rits_2int_10 <- gb_with_parts10 |>
  mutate(rit_output = map2(
    .x = nuc_id, .y = gbk,
    .f = ~{
      pb$tick()
      print(.x)
      rit_finder2(x = .x, gbk_xml = .y)
    }
  ))
beep()
rits_2int_10
rits_2int_10 |> unnest() |> unnest()
beep()
write_rds(rits_2int_10, './data/CDD/rit_finder_rs/rits_2int_10.rds')

rm(rits_2int_10, gb_with_parts10)


## 11 ----

gb_with_parts11 <- read_rds('./data/CDD/genbank/GB_w_parts_1001_1100.rds') |>
  unnest(gbk) |> 
  filter(!nuc_id %in% cause_error) |>
  mutate(size= map_dbl(gbk, ~object.size(.x))) |> 
  filter(size < 1.7e7)


pb <- make_pb(nrow(gb_with_parts11))
rits_2int_11 <- gb_with_parts11 |>
  mutate(rit_output = map2(
    .x = nuc_id, .y = gbk,
    .f = ~{
      pb$tick()
      print(.x)
      rit_finder2(x = .x, gbk_xml = .y)
    }
  ))
beep()
rits_2int_11
rits_2int_11 |> unnest() |> unnest()
beep()
write_rds(rits_2int_11, './data/CDD/rit_finder_rs/rits_2int_11.rds')
rm(rits_2int_11, gb_with_parts11)


## 12 -----

gb_with_parts12 <- read_rds('./data/CDD/genbank/GB_w_parts_1101_1200.rds') |>
  unnest(gbk) |> 
  filter(!nuc_id %in% cause_error) |>
  mutate(size= map_dbl(gbk, ~object.size(.x))) |> 
  filter(size < 1.7e7)


pb <- make_pb(nrow(gb_with_parts12))
rits_2int_12 <- gb_with_parts12 |>
  mutate(rit_output = map2(
    .x = nuc_id, .y = gbk,
    .f = ~{
      pb$tick()
      print(.x)
      rit_finder2(x = .x, gbk_xml = .y)
    }
  ))
beep()
rits_2int_12
rits_2int_12 |> unnest() |> unnest()
beep()
write_rds(rits_2int_12, './data/CDD/rit_finder_rs/rits_2int_12.rds')

rm(rits_2int_12, gb_with_parts12)


## 13 ------

gb_with_parts13 <- read_rds('./data/CDD/genbank/GB_w_parts_1201_1300.rds') |>
  unnest(gbk) |> 
  filter(!nuc_id %in% cause_error) |>
  mutate(size= map_dbl(gbk, ~object.size(.x))) |> 
  filter(size < 1.7e7)


pb <- make_pb(nrow(gb_with_parts13))
rits_2int_13 <- gb_with_parts13 |>
  mutate(rit_output = map2(
    .x = nuc_id, .y = gbk,
    .f = ~{
      pb$tick()
      print(.x)
      rit_finder2(x = .x, gbk_xml = .y)
    }
  ))
beep()
rits_2int_13
rits_2int_13 |> unnest() |> unnest()
beep()
write_rds(rits_2int_13, './data/CDD/rit_finder_rs/rits_2int_13.rds')

rm(rits_2int_13, gb_with_parts13)



## 14 -----

gb_with_parts14 <- read_rds('./data/CDD/genbank/GB_w_parts_1301_1400.rds') |>
  unnest(gbk) |> 
  filter(!nuc_id %in% cause_error) |>
  mutate(size= map_dbl(gbk, ~object.size(.x))) |> 
  filter(size < 1.7e7)

pb <- make_pb(nrow(gb_with_parts14))
rits_2int_14 <- gb_with_parts14 |>
  mutate(rit_output = map2(
    .x = nuc_id, .y = gbk,
    .f = ~{
      pb$tick()
      print(.x)
      rit_finder2(x = .x, gbk_xml = .y)
    }
  ))
beep()
rits_2int_14
rits_2int_14 |> unnest() |> unnest()
beep()
write_rds(rits_2int_14, './data/CDD/rit_finder_rs/rits_2int_14.rds')
rm(rits_2int_14, gb_with_parts14)

