# Start analyzing RITs that I found

# caluclate protein lengths and overlaps.
# verify output makes sense in most cases manually.
# extract sequences. reorient if on minus strand ('all reverse')
# see which are unique.
# which are duplicated within strains, among strains.
#

library(tidyverse)
library(janitor)
library(seqinr)

## Tidy up data =========================================
## and figure out which ids didn't work out...

# Results from rit_finder in P2_A5
rits <- Sys.glob('./data/CDD/rit_finder/*rs*.rds') |> 
  map(read_rds) |> 
  purrr::reduce(bind_rows) |> 
  unnest(rit_output, keep_empty = T)

# 98 of 517 missing CDS
missing_cds <- 
  rits |> filter(success == FALSE) |> 
  select(nuc_id, file) 
write_rds(missing_cds, './data/CDD/rit_finder/missing_cds.rds')

# files where genbank records came from
files <- rits |> 
  select(nuc_id, file)

# keep only searches that worked -> 419
# calculate integrase lengths
rit_data <- rits |> 
  anti_join(missing_cds, by = c("nuc_id", "file")) |> 
  select(-nuc_id_fail, -file) |> 
  unnest(cols = c(rits), keep_empty = T) |> 
  distinct()
glimpse(rit_data)

# 20 more nuc ids simply didn't return any results...
# just a tibble full of NAs - but 'success' flag not appropriate!
missing_results <- rit_data |> 
  filter(is.na(p1_id)) |> 
  select(nuc_id) |> 
  left_join(files, by = "nuc_id")
glimpse(missing_results)
write_rds(missing_results, './data/CDD/rit_finder/missing_results.rds')

# Results from nucleotides that worked; add unique rit_ids
rit_results <- rit_data |> 
  anti_join(missing_results,  by = "nuc_id") |> 
  select(-success, -contains('check')) |> 
  mutate(rit_id = paste0('RIT_', row_number(), '_', nuc_id)) |> 
  relocate(rit_id)

glimpse(rit_results)
rm(rit_data, rits)

## Add taxonomy ========================================

# create taxonomy table
rit_taxa <- read_rds('./data/CDD/id_data_fixed.rds') |> 
  select(nuc_id, tax_id) |> 
  distinct() |> 
  filter(nuc_id %in% rit_results$nuc_id) |> 
  left_join(read_rds('./data/CDD/tax_data_fixed.rds'),  by = "tax_id") |>
  unnest(tax_lineage) |> 
  clean_names()
glimpse(rit_taxa)

# count RITs by taxa, not sure about which elements are copied (not n unique but n total)
rit_results |> 
  left_join(rit_taxa) |> 
  dplyr::count(phylum, class, order, family, genus, 
               species, tax_id, nuc_id, name = 'n_tot') |> 
  arrange(desc(n_tot)) # |> View()

# get protein lengths; join taxonomy
rit_results <- rit_results |> 
  mutate(
    p1_length = (p1_stop - p1_start + 1) / 3,
    p2_length = (p2_stop - p2_start + 1) / 3,
    p3_length = (p3_stop - p3_start + 1)/ 3,
  ) |> 
  select(-contains('product')) |> 
  left_join(rit_taxa |> select(nuc_id, tax_id))

# 834 putative rit elements in 399 nucleotides; 187 unique taxa
skimr::skim(rit_results)


## Check properties ============================================

# orientations of RitA,B,C
rit_results |> 
  dplyr::count(rit_orientation)
# may the rit_finder is getting messed up by overlapping RIT elements... never frf or rfr, only fff, rrr, ffr, rff
rit_results |> 
  dplyr::count(p1_orientation, p2_orientation, p3_orientation)

# find RITs with not same orientation
non_consistent_orientation <- 
  rit_results |> 
  filter(rit_orientation == 'not same orientation')
# non_consistent_orientation |> View()

# 4/12 have oddly large overlaps, too
# -656 is gap of 656bp
# 98 is 
non_consistent_orientation |> 
  select(rit_id, contains('overlap'), contains('orientation'))

# find RITs with odd overlaps
odd_overlaps <- 
  rit_results |> 
  filter(overlap_p1p2 < -100 | overlap_p2p3 < -100) 
# odd_overlaps |> View()

# combine those RITs
weird_overlap_orientation <- 
  bind_rows(non_consistent_orientation, odd_overlaps)
rm(odd_overlaps, non_consistent_orientation)

# remove those RITs from main results... for inspection later 
# TODO inspection
rit_results <- rit_results |> 
  filter(!rit_id %in% weird_overlap_orientation$rit_id)

# 834 putative rit elements in 393 nucleotides; 187 unique taxa
skimr::skim(rit_results)

rit_results |> 
  dplyr::count(rit_orientation, p1_pred, p2_pred, p3_pred)

## Extract sequences =========================================
# get dna sequences for each putative element from nuc_data file

# fn takes a vector of ids and extract their seqs from file
extract_sequences <- function(file, ids){
  read_rds(file) |> 
    filter(nuc_id %in% ids)
}

# for all nuc_data_fixed files, get RIT sequences, take rev.comp if orientation of genes is 'all reverse'
rit_dna <- 
  Sys.glob('./data/CDD/nuc_data_fixed/*') |> 
  map(~extract_sequences(.x, ids = rit_results$nuc_id)) |> 
  purrr::reduce(bind_rows) |> 
  left_join(
    rit_results |> 
              transmute(rit_id, 
                        nuc_id, 
                        rit_start = p1_start,
                        rit_stop = p3_stop, 
                        rit_orientation),
    by = 'nuc_id'
  ) |> 
  mutate(
    rit_dna = str_sub(nuc_seq, rit_start, rit_stop),
    # get complement if RIT element on minus strand
    rit_dna = ifelse(
      test = rit_orientation == 'all reverse',
      yes = rit_dna |> 
        str_split('', simplify = T) |>
        seqinr::comp(ambiguous = T) |> 
        rev() |> 
        str_c(collapse = ''),
      no = rit_dna
    ),
    nuc_accession = str_trim(str_extract(nuc_name, '.*?\\s')),
    rit_length = rit_stop - rit_start + 1
  ) |> 
  relocate(rit_id, rit_length, rit_start, rit_stop)


## Time to check which sequences are unique...


rit_dna |> 
  select(-contains('seq')) |> 
  View()
  

