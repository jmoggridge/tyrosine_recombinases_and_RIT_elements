#### P2_A6 Start analyzing RITs that I found ====================

### Objectives ----

## Remove any nucleotide ids that didn't return any data
# - no genbank cds
# - just no rits

## Filter out any elements with large overlaps or gaps and those with mixed protein orientations - > element interrupted by another RIT element in several cases, not sure if unique...

## Extract RIT element DNA sequences. Reorient them if on minus strand (rit_orientation  = 'all reverse').
# verify output makes sense in most cases manually.

## EDA
# RIT lengths
# Plot protein lengths and overlaps.
# Flanking genes

## Taxonomy
# See how many & which are unique. Are copies shared within strains, among strains, range among taxa.

## Protein alignments  -> phylogeny of RIT elements / clustering
# *Translate protein sequences since we didn't take them from the genbank records....
# Identify neighboring genes & their distance


# Setup ------

library(tidyverse)
library(janitor)
library(seqinr)
library(glue)

# 07/20 Results from rit_finder
rits <- Sys.glob('./data/CDD/rit_finder_rs/*rs*.rds') |>
  map(read_rds) |>
  purrr::reduce(bind_rows) |>
  unnest(rit_output, keep_empty = T)

# all 521 nuc ids unique
length(unique(rits$nuc_id))

## Tidy up data ===============================================
## and figure out which ids didn't work out...

# 98 are missing CDS
missing_cds <- 
  rits |> filter(success == FALSE) |> 
  select(nuc_id, file) |> 
  mutate(ISSUE = 'missing CDS')
glimpse(missing_cds)
write_rds(missing_cds, './data/CDD/rit_finder/missing_cds.rds')

# files where genbank records came from
gb_files <- rits |> 
  select(nuc_id, file)

# keep only searches that worked -> 423 nuc ids w success = T
rit_data <- rits |> 
  anti_join(missing_cds, by = c("nuc_id", "file")) |> 
  select(-nuc_id_fail, -file)
glimpse(rit_data)

# there are 1076 elements when those 419 rows unnested
rit_data_unnest <- rit_data |> 
  unnest(cols = c(rits), keep_empty = T) |> 
  distinct()
glimpse(rit_data_unnest)


# These 22 nuc ids simply didn't return any results...
# just a tibble full of NAs - but 'success' flag not appropriate!
# some issue in parsing genbank data
missing_results <- rit_data_unnest |> 
  filter(is.na(p1_id)) |> 
  select(nuc_id) |> 
  left_join(gb_files, by = "nuc_id") |> 
  mutate(ISSUE = 'has CDS, rit_finder fail')
glimpse(missing_results)
write_rds(missing_results, './data/CDD/rit_finder/missing_results.rds')

# Results from nucleotides that worked -> 1054 possible elements
# add unique trio id to each row
rit_results <- rit_data_unnest |> 
  anti_join(missing_results,  by = "nuc_id") |> 
  select(-success, -contains('check')) |> 
  group_by(nuc_id) |> 
  mutate(trio_id = glue('trio_{row_number()}_{nuc_id}')) |> 
  ungroup() |> 
  relocate(trio_id)
glimpse(rit_results)

rm(rit_data, rit_data_unnest, rits, missing_cds, missing_results)

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

# count RITs by taxa, not sure about which elements are copied
# (not by unique elements, but by total elements)
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
    p3_length = (p3_stop - p3_start + 1) / 3,
  ) |> 
  select(-contains('product')) |> 
  left_join(rit_taxa |> select(nuc_id, tax_id))

# 834 putative rit elements in 399 nucleotides; 187 unique taxa
skimr::skim(rit_results)


## Check properties ============================================

# orientations of RitA,B,C; 101 not all same orientation
rit_results |> 
  dplyr::count(rit_orientation)

# may the rit_finder is getting messed up by overlapping RIT elements... never frf or rfr, only fff, rrr, ffr, rff
rit_results |> 
  dplyr::count(p1_orientation, p2_orientation, p3_orientation)

# find 101 RITs with not same orientation
non_consistent_orientation <- 
  rit_results |> 
  filter(rit_orientation == 'not same orientation')

# TODO which are overlapped by another RIT?
# a bunch of these are RIT ABC>/<BA or <CBA/BC> interrupted / duplicated?
overlapping_trios <-
  non_consistent_orientation |> 
  bind_rows(rit_results |> filter(nuc_id %in% non_consistent_orientation$nuc_id )) |> 
  arrange(nuc_id, trio_id, p1_start) |> 
  select(
    nuc_id, starts_with('trio_'), rit_orientation, ends_with('_pred'),
    p1_start, p1_stop, p2_start, p2_stop, p3_start, p3_stop,
    contains('overlap'), contains('orientation'),
    everything()
    ) |> 
  select(-c(contains('seq'), contains('flank'), contains('adjacent'),
            contains('cds'),
            matches('p[1-3]_id'))) |> 
  distinct() |> 
  group_by(nuc_id) |> 
  filter(lag(p3_stop) > p1_start | lead(p1_start) < p3_stop)

## Can see where there are 4+ integrases in a row in some cases
## quartets / quintets, often on different strands, not always though
## often contains a 'no consensus'-prediction for integrase. (degraded?)

# overlapping_trios |> View()

### Only a few with non-consistent orientation aren't overlapping with another trio

# non_consistent_orientation |> 
#   filter(!trio_id %in% overlapping_trios$trio_id) |> 
#   View()

# find 242 trios with odd overlaps between p1,p2,p3
odd_overlaps <- 
  rit_results |> 
  filter(overlap_p1p2 < -25 | overlap_p2p3 < -25 | 
           overlap_p1p2 > 25 | overlap_p2p3 > 25 ) |> 
  relocate(contains('overlap'))
odd_overlaps |> print()

# odd_overlaps |> View()

# combine those trios
weird_overlap_or_orientation <- 
  bind_rows(non_consistent_orientation, odd_overlaps)
weird_overlap_or_orientation
rm(odd_overlaps, non_consistent_orientation)


# remove those trios from main results... for inspection later 
# TODO inspection
# 899 likely RITs in trios
rit_candidates <- rit_results |> 
  filter(!trio_id %in% weird_overlap_or_orientation$trio_id) |> 
  group_by(nuc_id) |> 
  mutate(rit_id = glue('rit_{row_number()}_{nuc_id}')) |> 
  ungroup() |> 
  relocate(tax_id, nuc_id, rit_id) |> 
  select(-c(trio_id, contains('_adjacent'))) |> 
  arrange(tax_id, nuc_id, rit_id)

# rit_candidates |> View()

# 806 putative rit elements
skimr::skim(rit_candidates)

rit_candidates |> 
  dplyr::count(rit_orientation, p1_pred, p2_pred, p3_pred)

write_rds(weird_overlap_or_orientation,
          './data/CDD/rits_w_weird_overlaps_or_orientation.rds')
write_rds(overlapping_trios,
          './data/CDD/rit_overlapping_trios.rds')

rm(weird_overlap_or_orientation, overlapping_trios, gb_files)

## Extract sequences =========================================
# get dna sequences for each putative element from nuc_data file

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
  map(~{read_rds(.x) |> filter(nuc_id %in% rit_candidates$nuc_id)}) |> 
  purrr::reduce(bind_rows) |> 
  left_join(
    rit_candidates |> 
      transmute(
        tax_id, nuc_id,  rit_id, 
        rit_start = p1_start, rit_stop = p3_stop, 
        rit_orientation, p1_length, p2_length, p3_length,
        ),
    by = 'nuc_id'
  ) |> 
  mutate(
    
    strand = ifelse(rit_orientation == 'all reverse', 'minus', 'plus'),
    # extract rit sequences
    rit_dna = ifelse(
      test = strand == 'minus',
      yes = revcomp(str_sub(nuc_seq, rit_start, rit_stop)),
      no = str_sub(nuc_seq, rit_start, rit_stop)
      ),
    # extract upstream seq
    rit_dna_upstream = ifelse(
      test = strand == 'minus',
      yes = revcomp(str_sub(nuc_seq, rit_stop, rit_stop + 200)),
      no = str_sub(nuc_seq, rit_start - 200, rit_start)
    ),
    # extract downstream seq
    rit_dna_downstream = ifelse(
      test = strand == 'minus',
      yes = revcomp(str_sub(nuc_seq, rit_start - 200, rit_start)),
      no = str_sub(nuc_seq, rit_stop, rit_stop + 200)
    )
  ) |> 
  mutate(
    nuc_accession = str_trim(str_extract(nuc_name, '.*?\\s')),
    rit_length = rit_stop - rit_start + 1
  ) 

# Assign rit_dna_id to unique sequences
rit_candidates <- rit_dna |> 
  select(rit_dna) |> 
  distinct() |> 
  mutate(rit_dna_id = as_factor(row_number())) |> 
  right_join(rit_dna, by = "rit_dna") |> 
  left_join(rit_candidates, 
            by = c("nuc_id", "tax_id", "rit_id", "rit_orientation", 
                   "p1_length", "p2_length", "p3_length")) |> 
  select(-c(nuc_seq, element_CDS_length, matches('p[1-3]_orientation'))) |> 
  relocate(tax_id, nuc_id, rit_id, rit_dna_id,
           starts_with('rit_'), starts_with('nuc_'),
           starts_with('p1'), starts_with('p2'), starts_with('p3')
           ) |> 
  relocate(-ends_with('dna'), -ends_with('upstream'), -ends_with('downstream'))

glimpse(rit_candidates)

length(unique(rit_candidates$rit_dna))
rm(rit_dna)

rit_candidates |> 
  relocate(contains('pred')) |> 
  View()


## triaging obs for inspection...

# pairs of nucleotide records from same taxa, same regions, have same trios
identical_nucleotides <- 
  c('1597162718', # 1601754439
    ''
    )

# element that overlaps a RIT, where flanking gene is also an integrase
has_integrase_beside_RIT <- rit_candidates |> 
  filter(rit_id %in% c(
    'rit_6_1680973734', 'rit_1_1308917638', 'rit_3_1597162718',
    'rit_4_1597162718'
    ))

# seem to be legit RIT elements (filter out to simplify set to review)
checked_rits_good <-  rit_candidates |> 
  filter(rit_id %in% c(
    'rit_1_218292240', 'rit_1_806916613', 'rit_2_806916613',
    'rit_1_1545075283', 'rit_2_1545075283', 'rit_1_1258461527', 'rit_1_1253857670',
    'rit_1_47118328', 'rit_2_47118328', 'rit_3_47118328', 
    'rit_1_1207197702', 'rit_1_1198444759'
    
    ))

rit_candidates |> 
  filter(!rit_id %in% has_integrase_beside_RIT$rit_id) |> 
  filter(!rit_id %in% checked_rits_good$rit_id) |> 
  filter(!nuc_id %in% identical_nucleotides) |> 
  relocate(contains('pred')) |> 
  select(p1_pred:overlap_p2p3) |> 
  View()




## Check which sequences are unique...

## CHECKING RIT ELEMENT DNA SEQUENCES FOR ACCURACY - MATCHING NCBI data


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

## EDA ------------------------------------------------------------

## Plot frequency histogram of unique RIT sequence
# 420 unique sequences; some as many as 16
rit_candidates |> 
  dplyr::count(rit_dna_id) |> arrange(desc(n)) |> 
  ggplot(aes(n, fill = rit_dna_id)) +
  geom_histogram(position = 'stack', color = 'darkgray', size = 0.15,
                 show.legend = F) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  labs(title = 'Frequency of unique RIT elements in dataset of 806 sequences',
       x = 'number of unique RITs elements') +
  theme_classic() 



# get lengths of RitA,B,C
rit_candidates <- rit_dna |> 
  relocate(contains('rit'), contains('nuc')) |> 
  mutate(
    strand = ifelse(rit_orientation == 'all forward', '+', '-'),
    rit_A_length = ifelse(strand == '+', p1_length, p3_length),
    rit_B_length = p2_length,
    rit_C_length = ifelse(strand == '+', p3_length, p1_length)
  ) |> 
  select(-c(p1_length, p2_length, p3_length))

rit_candidates |> 
  ggplot(aes(nchar(rit_dna))) +
  geom_histogram() +
  labs(x = 'length (bp)', 
       title = 'RIT element size distribution')
# 
# rit_candidates |> 
#   pivot_longer(cols = c(rit_A_length, rit_B_length, rit_C_length)) |> 
#   mutate(name = str_remove(name, '_length') |> str_replace('r', 'R')) |> 
#   ggplot(aes(x = value)) +
#   geom_histogram(show.legend = F) +
#   facet_grid(.~name) +
#   labs(x = 'length (AA)')


