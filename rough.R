
library(tidyverse)
library(janitor)
library(seqinr)
library(glue)


## Functions ============================================================

# gets reverse complement of dna sequence
revcomp <- function(dna){
  map_chr(dna, ~{
    .x |> str_split('', simplify = T) |>
      seqinr::comp(ambiguous = T, forceToLower = F) |> 
      rev() |> 
      str_c(collapse = '') })
}

## Tidy up data =========================================================
## and figure out which ids didn't work out...

# Dataset - 07/20 Results from rit_finder
rits <- Sys.glob('./data/CDD/rit_finder_rs/*rs*.rds') |>
  map(read_rds) |>
  purrr::reduce(bind_rows) |>
  unnest(rit_output, keep_empty = T)
glimpse(rits)
# all 521 nuc ids unique
length(unique(rits$nuc_id))

# 98 are missing CDS
missing_cds <- 
  rits |> filter(success == FALSE) |> 
  select(nuc_id, file) |> 
  mutate(ISSUE = 'missing CDS')
glimpse(missing_cds)
write_rds(missing_cds, './data/CDD/rit_finder/missing_cds.rds')

# keep only searches that worked -> 423 nuc ids w success = T
# unnesting: there are 1076 integrase trios in those 419 nucleotide records 
rits <- rits |> 
  anti_join(missing_cds, by = c("nuc_id", "file")) |> 
  unnest(cols = c(rits), keep_empty = T)
glimpse(rits)

# files where genbank records came from
gb_files <- rits |> select(nuc_id, file)

# These 22 nuc ids simply didn't return any results...
# just a tibble full of NAs - but 'success' flag not appropriate!
# some issue in parsing genbank data
missing_results <- rits |> 
  filter(is.na(p1_id) | is.na(p2_id) | is.na(p3_id)) |> 
  select(nuc_id) |> 
  left_join(gb_files, by = "nuc_id") |> 
  mutate(ISSUE = 'has CDS, rit_finder fail')
glimpse(missing_results)
write_rds(missing_results, './data/CDD/rit_finder/missing_results.rds')

# Results from nucleotides that worked -> 1054 possible elements
# add unique trio id to each row
rit_rs <- rits |> 
  anti_join(missing_results,  by = "nuc_id") |> 
  select(-c(file, success, nuc_id_fail, element_CDS_length, 
            rit_dist_check, rit_overlap_check, rit_length_check,
            rit_integrase_check, contains('adjacent'))) |> 
  group_by(nuc_id) |> 
  mutate(trio_id = glue('trio_{row_number()}_{nuc_id}')) |> 
  ungroup() |> 
  relocate(trio_id) 
glimpse(rit_rs)

rm(rits, missing_cds, missing_results, gb_files)

## Link to taxonomy =======================================================

# create taxonomy table for 401 associated tax_ids
rit_taxa <- 
  read_rds('./data/CDD/id_data_fixed.rds') |> 
  select(nuc_id, tax_id) |> 
  distinct() |> 
  filter(nuc_id %in% rit_rs$nuc_id) |> 
  left_join(read_rds('./data/CDD/tax_data_fixed.rds'),  by = "tax_id") |>
  unnest(tax_lineage) |> 
  clean_names()
glimpse(rit_taxa)

# count trios by taxa, not sure about which elements are copied yet
# (not by unique elements, but by total elements)
rit_rs |> 
  left_join(rit_taxa) |> 
  dplyr::count(phylum, class, order, family, genus, 
               species, tax_id, nuc_id, name = 'n_tot') |> 
  arrange(desc(n_tot)) # |> View()

# link a taxonomy id to each nucleotide id
rit_rs <- rit_rs |> 
  left_join(rit_taxa |> select(nuc_id, tax_id)) |> 
  relocate(tax_id, nuc_id, trio_id)
glimpse(rit_rs)



## Process sequences ========================================================
# get dna sequences for each putative element from nuc_data file
# extract trio sequence, upstream and downstream 200bp

# for all nuc_data_fixed files, read in genbank sequence and join to rit results.
rit_dna <- 
  Sys.glob('./data/CDD/nuc_data_fixed/*') |> 
  map(~{read_rds(.x) |> filter(nuc_id %in% rit_rs$nuc_id)}) |> 
  purrr::reduce(bind_rows) |> 
  left_join(rit_rs |> 
              transmute(
                tax_id, nuc_id,  trio_id, 
                rit_start = p1_start, 
                rit_stop = p3_stop, 
                rit_orientation,
              ),
            by = 'nuc_id'
  )

glimpse(rit_dna)

# Generate sequence features
rit_dna <- 
  rit_dna |> 
  mutate(
    nuc_length = nchar(nuc_seq),
    # trio length, strand, parent record.
    rit_length = rit_stop - rit_start + 1,
    rit_strand = ifelse(rit_orientation == 'all reverse', 'minus', 'plus'),
    nuc_accession = str_trim(str_extract(nuc_name, '.*?\\s')),
    # extract trio dna sequence
    rit_dna = ifelse(
      test = rit_strand == 'minus',
      yes = revcomp(str_sub(nuc_seq, rit_start, rit_stop)),
      no = str_sub(nuc_seq, rit_start, rit_stop)
    )) |> 
  rowwise() |> 
  mutate(
    # can't try to substring past (1:length) range
    pl_upstream_bgn = max(1, rit_start - 200),
    pl_downstream_end = min(nuc_length, rit_stop + 200),
    min_upstream_bgn = min(nuc_length, rit_stop + 200),
    min_downstream_end = max(1, rit_stop + 200),
    # extract upstream seq
    rit_dna_upstream = ifelse(
      test = rit_strand == 'minus',
      yes = revcomp(str_sub(nuc_seq, rit_stop, min_upstream_bgn)),
      no = str_sub(nuc_seq, pl_upstream_bgn, rit_start)
    ),
    # extract downstream seq
    rit_dna_downstream = ifelse(
      test = rit_strand == 'minus',
      yes = revcomp(str_sub(nuc_seq, min_downstream_end, rit_start)),
      no = str_sub(nuc_seq, rit_stop, pl_downstream_end)
    )
  ) |> 
  ungroup() |> 
  mutate(across(.cols = matches('dna_(upstream|downstream)'),
                .fns = ~ifelse(nchar(.x)==1, NA, .x)))

glimpse(rit_dna) 
rit_dna |> 
  select(-matches('_bgn|_end|_seq|nuc_length')) |> 
  select(rit_orientation, matches('pl|mn|start|stop|dna_up|dna_down')) |> 
  View()

# join sequences data to main results data
# add protein lengths in AAs from nucleotide ranges
rit_rs <- rit_rs |> 
  left_join(rit_dna) |> 
  mutate(
    p1_length = (p1_stop - p1_start + 1) / 3,
    p2_length = (p2_stop - p2_start + 1) / 3,
    p3_length = (p3_stop - p3_start + 1) / 3
    ) |> 
  select(-nuc_seq)

## flip p1 & p3 for trios on the minus strand.
rit_rs_aligned <- rit_rs |>
  filter(rit_strand == 'minus') |> 
  rename_with(~str_replace(.x, 'p1_', 'r3_')) |> 
  rename_with(~str_replace(.x, 'p3_', 'r1_')) |> 
  rename_with(~str_replace(.x, 'r1_', 'p1_')) |>
  rename_with(~str_replace(.x, 'r3_', 'p3_')) |>
  rename_with(~str_replace(.x, 'lflank_.', 'xflk_')) |> 
  rename_with(~str_replace(.x, 'rflank_', 'yflk_')) |> 
  rename_with(~str_replace(.x, 'xflk_', 'rflank_')) |> 
  rename_with(~str_replace(.x, 'yflk_', 'lflank_')) |> 
  relocate(tax_id, nuc_id, trio_id, contains('p1'), contains('p2'), contains('p3')) |> 
  bind_rows(rit_rs |> filter(rit_strand == 'plus'))

rit_rs |> View()
rit_rs_aligned |> View()
skimr::skim(rit_rs)
rm(rit_dna)



## FILTERING TRIOS ======================================================
## Filter trios (1054) to retain only RIT elements:
# 1. remove incorrect orientations

# orientations of RitA,B,C; 101 not all same orientation
rit_rs |> dplyr::count(rit_orientation)
# may the rit_finder is getting messed up by overlapping RIT elements... never frf or rfr, only fff, rrr, ffr, rff
rit_rs |> dplyr::count(p1_orientation, p2_orientation, p3_orientation)

# from all integrase trios (1054) keep only trios with three integrases
# with consistent orientation (953)
oriented_trios <- rit_rs |> 
  filter(rit_orientation %in% c('all forward', 'all reverse'))
glimpse(oriented_trios)

# rit_elements <- oriented_trios |> 
  
# Those with non-consistent orientation are
# mostly partial elements from strains with complete elements (101)
# there are 71 unique nuc ids, but 70 are also in rit_elements!
non_oriented <- rit_rs |> filter(!trio_id %in% rit_elements$trio_id)
non_oriented
length(non_oriented |> pull(nuc_id) |> unique())
sum(non_oriented |> pull(nuc_id) |> unique() %in% rit_elements$nuc_id)

# non_oriented ones
View(
  non_oriented |> 
    select(-contains('cds'), -contains('dna'), rit_strand) |> 
    relocate(tax_id, nuc_id, trio_id, 
             matches('p[1-3]_pred'), matches('p1_st'),  matches('p2_st'), 
             matches('p3_st'), contains('rit'))
)

# see where these are next to a complete element

# Examine elements that don't have the RIT A-B-C arrangement (131/953)
non_ABC <- rit_elements |> 
  filter(!(p1_pred == 'RitA' & p2_pred == 'RitB' & p3_pred == 'RitC')) |> 
  filter(!(p1_pred == 'RitC' & p2_pred == 'RitB' & p3_pred == 'RitA')) |> 
  select(-contains('cds'), -contains('dna'), rit_strand) |> 
  relocate(tax_id, nuc_id, trio_id, 
           matches('p[1-3]_pred'), matches('p1_st'),  matches('p2_st'), 
           matches('p3_st'), contains('rit'))

# any that are ACB or BCA? 20
ACB_BCA <- non_ABC |> 
  rowwise() |> 
  filter(
    all(p1_pred == 'RitA', p2_pred == 'RitC', p3_pred == 'RitB') |
    all(p1_pred == 'RitB', p2_pred == 'RitC', p3_pred == 'RitA')
    )
ACB_BCA

# Any that are CAB or BAC? 
CAB_BAC <- non_ABC |> 
  filter(
    all(p1_pred == 'RitC', p2_pred == 'RitA', p3_pred == 'RitB') |
    all(p1_pred == 'RitB', p2_pred == 'RitA', p3_pred == 'RitC')
  )
CAB_BAC  

# which have flanking Rit A or C protein & a non-Rit at p1 or p3?


non_ABC |> View()
  
rit_elements |> 
  select(-contains('cds'), -contains('dna'), rit_strand) |> 
  relocate(tax_id, nuc_id, trio_id, 
           matches('p[1-3]_pred'), matches('p1_st'),  matches('p2_st'), 
           matches('p3_st'), contains('rit')) |> 
  View()
  