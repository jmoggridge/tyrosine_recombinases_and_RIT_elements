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

# 806 putative rit elements
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

# TODO unvectorize!!
# gets reverse complement of dna sequence
revcomp <- function(dna){
  map_chr(dna, ~{
    .x |> str_split('', simplify = T) |>
    seqinr::comp(ambiguous = T, forceToLower = F) |> 
    rev() |> 
    str_c(collapse = '') })
}

# for all nuc_data_fixed files, read in genbank sequence and join to rit results.
rit_dna <- 
  Sys.glob('./data/CDD/nuc_data_fixed/*') |> 
  map(~extract_sequences(.x, ids = rit_results$nuc_id)) |> 
  purrr::reduce(bind_rows) |> 
  left_join(
    rit_results |> 
      transmute(
        rit_id, 
        nuc_id, 
        rit_start = p1_start,
        rit_stop = p3_stop, 
        rit_orientation,
        p1_length, p2_length, p3_length,
        ),
    by = 'nuc_id'
  ) |> 
  # extract rit sequences
  mutate(rit_dna_in_nuc = str_sub(nuc_seq, rit_start, rit_stop)) |> 
  # get complement if RIT element on minus strand
  mutate(
    rit_dna = ifelse(
      test = rit_orientation == 'all reverse',
      yes = revcomp(rit_dna_in_nuc),
      no = rit_dna_in_nuc
    ),
    nuc_accession = str_trim(str_extract(nuc_name, '.*?\\s')),
    rit_length = rit_stop - rit_start + 1
  ) |> 
  relocate(rit_id, contains('rit'), contains('nuc'), contains('dna'))


rit_dna <- rit_dna |> 
  select(rit_dna) |> 
  distinct() |> 
  mutate(rit_dna_id = as_factor(row_number())) |> 
  right_join(rit_dna)

# rit_dna |> View()

# 421 unique sequences; some as many as 16
rit_dna_df2 |> 
  dplyr::count(rit_dna_id) |> arrange(desc(n)) |> 
  ggplot(aes(n, fill = rit_dna_id)) +
  geom_histogram(position = 'stack', color = 'white', size = 0.05,
                 show.legend = F) +
  labs(title = 'Frequency of unique RIT elements in dataset of 806 sequences',
       x = 'number of unique RITs elements') +
  theme_classic() 

length(unique(rit_dna$rit_dna))
# 
# unique_rit_df <- rit_dna |> 
#   left_join(unique_rit_dna)
#   arrange(rit_dna) |> 
#   group_by(rit_dna) |> 
#   nest(occurences = c(rit_id, rit_start, rit_stop, rit_orientation, nuc_id, nuc_name, nuc_accession)) |> 
#   ungroup() |> 
#   mutate(rit_dna_id = as.numeric(as_factor(rit_dna))) 
# unique_rit_df |> 
#   View()



# get lengths of RitA,B,C
rit_rs <- rit_dna |> 
  relocate(contains('rit'), contains('nuc')) |> 
  mutate(
    strand = ifelse(rit_orientation == 'all forward', '+', '-'),
    rit_A_length = ifelse(strand == '+', p1_length, p3_length),
    rit_B_length = p2_length,
    rit_C_length = ifelse(strand == '+', p3_length, p1_length)
  ) |> 
  select(-c(p1_length, p2_length, p3_length))

rit_rs |> 
  ggplot(aes(nchar(rit_dna))) +
  geom_histogram()

rit_rs |> 
  pivot_longer(cols = c(rit_A_length, rit_B_length, rit_C_length)) |> 
  mutate(name = str_remove(name, '_length') |> str_replace('r', 'R')) |> 
  ggplot(aes(x = value)) +
  geom_histogram(show.legend = F) +
  facet_grid(.~name) +
  labs(x = 'length (AA)')
## Time to check which sequences are unique...
