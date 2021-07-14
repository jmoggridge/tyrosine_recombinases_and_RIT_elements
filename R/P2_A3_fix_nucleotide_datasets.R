## P2_A3_fix_nucleotide_datasets
# deal with sections where n(id) != n(seq)
# rewrite files as unnested data with _fixed.rds suffix
library(tidyverse)
library(glue)
library(beepr)
library(rentrez)
source('./R/P2_entrez_functions.R')

## Fix nuc data function -----

## tidy_up_nuc_set(nuc_file):

# open file, split good and bad sets -
# unnest the good ones,  remove any ids that are no longer 
# repair bad ones by joining nuc_id from summary using 'caption' - the accession code
# bind rows and resave _fixed.rds
# return dataframe of ids to recreate the index...

tidy_up_nuc_set <- function(nuc_file, nuc_summary, out_folder) {
  message(glue('Working on {nuc_file}'))
  message('original data:')
  
  # read dataset
  nuc_dat <- read_rds(nuc_file) |> 
    mutate(n_id = map_int(id, length),
           n_name = map_int(nuc_name, length),
           n_seq = map_int(nuc_seq, length)
    )
  print(nuc_dat)
  # unnest tidy data where possible
  nuc_tidy <- 
    nuc_dat |> 
    filter(n_id == n_name & n_id == n_seq) |> 
    transmute(nuc_id = id, nuc_name, nuc_seq) |> 
    unnest(cols = everything()) |>
    # discard any ids already filtered in P2
    filter(nuc_id %in% id_data$nuc_id)
  
  # deal with sections where n(id) != n(seq)
  # drop id, unnest the rest, extract 'caption' to match nuc_summary
  # inner_join seqs with nuc_summary to remove any data that was filtered out 
  # in P2_A2
  nuc_cleaned <- 
    nuc_dat |> 
    filter(n_id != n_name | n_id != n_seq) |> 
    transmute(nuc_id = id, nuc_name, nuc_seq) |> 
    select(-nuc_id) |> 
    unnest(everything()) |> 
    mutate(caption = str_extract(nuc_name, '^.*?\\.') |> 
             str_remove('\\.')) |> 
    relocate(caption) |> 
    inner_join(nuc_summary |> select(caption, nuc_id), by = 'caption') |> 
    select(nuc_id, nuc_name, nuc_seq)
  
  # recombine and write _filtered
  fixed_data <- 
    bind_rows(nuc_tidy, nuc_cleaned)
  message(glue('Fixed data has {nrow(fixed_data)} rows.\n\n'))
  
  new_file <- nuc_file |> 
    str_replace('./data/CDD/nuc_data/', out_folder)
  write_rds(fixed_data, file = new_file)
  
  return(fixed_data |> select(nuc_id))
}

## Main ----

id_data <- read_rds('./data/CDD/id_data_fixed.rds')
nuc_summary <- read_rds('./data/CDD/nuc_summary_fixed.rds')

in_folder <- './data/CDD/nuc_data/'
out_folder <- './data/CDD/nuc_data_fixed/'
dir.create(out_folder)

nuc_id_sets <- 
  Sys.glob(glue('{in_folder}*.rds')) |> 
  map(~tidy_up_nuc_set(nuc_file = .x, nuc_summary, out_folder))

# all the ids that I have records for already
nuc_ids_downloaded <- reduce(nuc_id_sets, bind_rows)

# nuc_ids missing data
missing_data <- id_data |> 
  filter(!nuc_id %in% nuc_ids_downloaded$nuc_id) 

# get any missing records for ids in id_data
new_records <- 
  missing_data |> 
  select(nuc_id) |>   
  mutate(
    token = map(
      nuc_id, 
      ~entrez_search(db = 'nuccore', term = .x, use_history = T)
    ),
    fasta = map(
      token,
      ~entrez_fetch(
        db = 'nuccore', 
        web_history = .x$web_history, 
        rettype = 'fasta'
      ))
  ) |> 
  mutate(ss = map(fasta, fasta_to_DNAss),
         nuc_name = map(ss, names),
         nuc_seq = map(ss, paste)) |> 
  select(-ss)

new_records |> 
  select(-token) |> unnest()
