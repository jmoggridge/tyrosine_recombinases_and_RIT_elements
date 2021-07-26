## P2_A6 -
# Clean up results from rit_finder table in P2_A5.

# join dna sequences, + upstream & downstream 200bp to 
#

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
## and results for certain ids didn't work out...
# These are labelled with success=FAIL tag.

# Some others just have no data, which could be a parsing failure somewhere or more likely that there is just no RIT element in the nucleotide.

# Dataset - 07/20 (.._rs_.. files) & 25 (rits_2...) Results from rit_finder
# 1887 obs of probable RIT elements
rits <- 
  Sys.glob('./data/CDD/rit_finder_rs/*.rds') |>
  map(read_rds) |>
  purrr::reduce(bind_rows) |>
  unnest(rit_output, keep_empty = T) |> 
  select(-gbk, -size, -tax_id)

# 1887 rits
# 1842 nuc ids 
glimpse(rits)
length(unique(rits$nuc_id))

# 98 rows are missing CDS --- but I added the replacements where possible (27)
missing_cds <- 
  rits |> filter(success == FALSE) |> 
  select(nuc_id, file) |> 
  mutate(ISSUE = 'missing CDS')
glimpse(missing_cds)

# keep only searches that worked -> 423 nuc ids w success = T
# unnesting: there are 1076 integrase trios in those 419 nucleotide records 
rits <- rits |> 
  anti_join(missing_cds, by = c("nuc_id", "file")) |> 
  unnest(cols = c(rits), keep_empty = T)
glimpse(rits)

# remove the ones that have been replaced from the missing list.
missing_cds <- missing_cds |> anti_join(rits, by = 'nuc_id')
write_rds(missing_cds, './data/CDD/missing_cds.rds')

# managed to fix 27...
# sum(missing_cds$nuc_id %in% rits$nuc_id)

# files where genbank records came from
gb_files <- rits |> select(nuc_id, file)

# These 166 nuc ids simply didn't return any results...
# just a tibble full of NAs - but 'success' flag not appropriate!
# some issue in parsing genbank data
missing_results <- rits |> 
  filter(is.na(p1_id) | is.na(p2_id) | is.na(p3_id)) |> 
  select(nuc_id) |> 
  left_join(gb_files, by = "nuc_id") |> 
  mutate(ISSUE = 'has CDS, rit_finder fail')
glimpse(missing_results)
write_rds(missing_results, './data/CDD/missing_results.rds')

# Results from nucleotides that worked -> 2596 possible elements
# add unique 'trio_id' to each row to track
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

# ## Join missing results ====================================
# 
# fixed_cds <- read_rds('data/CDD/RIT_gbk_noCDS_fixed.rds')
# glimpse(fixed_cds)

## Taxon links =======================================================

# create taxonomy table for 1623 associated tax_ids
rit_taxa <- 
  read_rds('./data/CDD/id_data_fixed.rds') |> 
  select(nuc_id, tax_id) |> 
  distinct() |> 
  filter(nuc_id %in% rit_rs$nuc_id) |> 
  left_join(read_rds('./data/CDD/tax_data_fixed.rds'),  by = "tax_id") |>
  unnest(tax_lineage) |> 
  clean_names() |> 
  select(-no_rank)
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



## DNA sequences ========================================================

# get dna sequences for each putative element from nuc_data file
# extract trio sequence, upstream and downstream 200bp
# flip p1/p3 when element is on minus strand so that p1 is first in order of transcription

##

# for all nuc_data_fixed files, read in genbank sequence and join to rit results.
rit_dna <- 
  Sys.glob('./data/CDD/nuc_data_fixed/*') |> 
  map(~{read_rds(.x) |> filter(nuc_id %in% rit_rs$nuc_id)}) |> 
  purrr::reduce(bind_rows) |> 
  left_join(
    rit_rs |> 
      transmute(
        tax_id, nuc_id,  trio_id, 
        rit_start = p1_start, 
        rit_stop = p3_stop, 
        rit_orientation),
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
    min_downstream_end = max(1, rit_start - 200),
    # extract upstream seq
    rit_dna_upstream = ifelse(
      test = rit_strand == 'minus',
      yes = revcomp(str_sub(nuc_seq, rit_stop+1, min_upstream_bgn)),
      no = str_sub(nuc_seq, pl_upstream_bgn, rit_start-1)
    ),
    # extract downstream seq
    rit_dna_downstream = ifelse(
      test = rit_strand == 'minus',
      yes = revcomp(str_sub(nuc_seq, min_downstream_end, rit_start-1)),
      no = str_sub(nuc_seq, rit_stop+1, pl_downstream_end)
    )
  ) |> 
  ungroup() |> 
  mutate(across(.cols = matches('rit_dna_up|rit_dna_down'),
                .fns = ~ifelse(nchar(.x)==1, NA, .x))) |> 
  select(-matches('_bgn|_end|_seq|nuc_length')) 

glimpse(rit_dna) 


# join sequences data to main results data
# add protein lengths in AAs from nucleotide ranges
rit_rs <- rit_rs |> 
  left_join(rit_dna) |> 
  mutate(
    p1_length = (p1_stop - p1_start + 1) / 3,
    p2_length = (p2_stop - p2_start + 1) / 3,
    p3_length = (p3_stop - p3_start + 1) / 3
    ) 

## flip p1 & p3 columns for trios on the minus strand so that p1 is the first
## gene in order of transcription
rit_rs_aligned <- rit_rs |>
  filter(rit_strand == 'minus') |> 
  rename_with(~str_replace(.x, 'p1_', 'r3_')) |> 
  rename_with(~str_replace(.x, 'p3_', 'r1_')) |> 
  rename_with(~str_replace(.x, 'r1_', 'p1_')) |>
  rename_with(~str_replace(.x, 'r3_', 'p3_')) |>
  rename_with(~str_replace(.x, 'lflank_', 'xflk_')) |> 
  rename_with(~str_replace(.x, 'rflank_', 'yflk_')) |> 
  rename_with(~str_replace(.x, 'xflk_', 'rflank_')) |> 
  rename_with(~str_replace(.x, 'yflk_', 'lflank_')) |> 
  relocate(tax_id, nuc_id, trio_id, contains('p1'), contains('p2'), contains('p3')) |> 
  bind_rows(rit_rs |> filter(rit_strand == 'plus'))

skimr::skim(rit_rs_aligned)
rm(rit_dna, rit_rs)

## Proteins set ================================================

# Create a collection of all the integrases retrieved from any trios
proteins_set <- rit_rs_aligned |> 
  pivot_longer(matches('p[1-3]_|[rl]flank_'), 
               names_to = c("set", ".value"),
               names_pattern = "(.+)_(.+)"
  ) |> 
  select(nuc_id, tax_id, trio_id, id, pred, product, 
         start, stop, orientation, length) |> 
  arrange(tax_id, trio_id, start) |> 
  rename(prot_id = id)
  
glimpse(proteins_set)
write_rds(proteins_set, 'data/CDD/rits_proteins_set.rds')
rm(proteins_set)

## FILTERING TRIOS ======================================================
## Filter trios to retain only RIT elements:
# 1. remove incorrect orientations

# orientations of RitA,B,C; some have not all same orientation
rit_rs_aligned |> dplyr::count(rit_orientation)
# may the rit_finder is getting messed up by overlapping RIT elements... never frf or rfr, only fff, rrr, ffr, rff
rit_rs_aligned |> dplyr::count(p1_orientation, p2_orientation, p3_orientation)

# from all integrase trios, keep only trios with three integrases
# with consistent orientation
oriented_trios <- rit_rs_aligned |> 
  filter(rit_orientation %in% c('all forward', 'all reverse'))
glimpse(oriented_trios)

# View(oriented_trios |> select(-matches('dna')) |> relocate(contains('_pred')))

# Most have  A-B-C arrangement
rit_abc_elements <- oriented_trios |>
  filter(p1_pred == 'RitA' & p2_pred == 'RitB' & p3_pred == 'RitC')
rit_abc_elements

# CTnDOT proteins upstream of a rit element in two nuc from a strain,
# annotated slightly differently
CTn_DOT_elements <-  oriented_trios |>
  rowwise() |> 
  filter(sum(c(p1_pred,p2_pred,p3_pred) == 'Int_CTnDOT') ==2)
CTn_DOT_elements

# some 'XIT element' with overlaps!
xer_elements <-  oriented_trios |>
  rowwise() |> 
  filter(sum(c(p1_pred,p2_pred,p3_pred) == 'Xer') >1)
xer_elements

# other trios with no protein already present in ABC elements
# mostly degraded RIT elements, have a lot of no-consensus predictions...
# these are probably partial domains with low scores across many subfamilies
# have unusual overlaps, suggestion degradation too.
rit_abc_proteins <- c(rit_abc_elements |> select(p1_id, p3_id) |> 
                        unlist())
non_overlapping_trios <- oriented_trios |>
  filter(!trio_id %in% rit_abc_elements$trio_id) |>   
  filter(!trio_id %in% xer_elements$trio_id) |> 
  filter(!p1_id %in% rit_abc_proteins & 
           !p2_id %in% rit_abc_proteins &
           !p3_id %in%rit_abc_proteins)
non_overlapping_trios
non_overlapping_trios |> dplyr::count(nuc_id) |> nrow()

# non_overlapping_trios |> relocate(lflank_pred, ends_with('_pred'), rflank_pred) |> View()

## Any non-overlapping that aren't degraded?? no, mostly degraded elements
non_overlapping_trios |> 
  filter(!trio_id %in% rit_abc_elements$trio_id) |>   
  select(matches('pred|id')) |> 
  rowwise() |> 
  filter(!any('no consensus' %in% c(p1_pred, p2_pred, p3_pred))) 
# Neither has RitABC
# Another is just degraded
# trio_2_1719553555 has a gap of 429 bp between RitA - RitB

## overlapping rit_abc element -> remove from results and ignore
overlapping_abc_trios <- oriented_trios |>
  filter(!trio_id %in% rit_abc_elements$trio_id) |>   
  filter(!trio_id %in% xer_elements$trio_id) |> 
  filter(p1_id %in% rit_abc_proteins | 
           p2_id %in% rit_abc_proteins | 
           p3_id %in% rit_abc_proteins)
overlapping_abc_trios
# overlapping_abc_trios |> View()


### account for triaged integrase trios ====================================

# 1260 integrase trios with consistent orientation
nrow(oriented_trios)
# 1113 good rit elements
nrow(rit_abc_elements)
# 1 weird Xer trio, to exclude from RITs results set
nrow(xer_elements)
# 146 left
rem <- nrow(oriented_trios) - nrow(rit_abc_elements) - nrow(xer_elements)
rem
# 121 are trios that overlap a rit ABC element
nrow(overlapping_abc_trios)
# 25 are not overlapping rit ABC elements - mostly degraded
nrow(non_overlapping_trios)
# All accounted for!
rem - nrow(non_overlapping_trios) - nrow(overlapping_abc_trios)


### Combine RIT elements, tidy table up
rit_elements <- 
  bind_rows(rit_abc_elements) |> 
  mutate(rit_id = row_number()) |> 
  mutate(p1_overlap = overlap_p1p2,
         p2_overlap = overlap_p2p3) |> 
  select(-starts_with('overlap'), -contains('check'), -rit_orientation) |> 
  pivot_longer(matches('p[1-3]|[rl]flank'), 
               names_to = c("prot_pos", ".value"),
               names_pattern = "(.+)_(.+)") |> 
  rename_with(.cols = c(id, pred, product, start, stop, 
                        orientation, length, overlap),
              .fn = ~paste0('prot_', .x)) |> 
  nest(protein_df = starts_with('prot_')) |> 
  relocate(ends_with('id'), starts_with('rit'))

glimpse(rit_elements)
write_rds(rit_elements, './results/rit_elements.rds')


# tidy up
rm(CTn_DOT_elements, rit_rs_aligned, overlapping_abc_trios, rit_abc_elements,
   non_overlapping_trios, oriented_trios, xer_elements,
   rit_abc_proteins, rem)



## Which are distinct elements? ==============================================

# add a rit_unique_seq_id to distinct sequences
rit_elements <- rit_elements |> 
  select(rit_dna) |> 
  distinct() |> 
  mutate(rit_unique_seq_id = glue('RIT_{row_number()}')) |> 
  right_join(rit_elements, by = 'rit_dna')

rit_elements |>
  dplyr::count(rit_unique_seq_id) |> 
  ggplot(aes(n)) +
  geom_histogram() +
  labs(x = 'Exact copies found in results', y = 'n distinct RIT elements',
       title = 'Distribution of the number of identical RIT sequences across all records')

# distinct elements by tax id, RIT sequence... 485.
distinct_rits <- rit_elements |> 
  group_by(tax_id, rit_unique_seq_id) |> 
  nest(rit_occurrences = c(
    nuc_id, nuc_accession, nuc_name, rit_id, trio_id,
    rit_length, protein_df, matches('rit_st|upstream|downstream'))
  ) |> 
  ungroup() |> 
  mutate(copies = map_int(rit_occurrences, ~nrow(.x))) |> 
  left_join(rit_taxa |> select(-nuc_id) |> distinct(), by = 'tax_id')
  
glimpse(distinct_rits)

write_rds(distinct_rits, 'results/distinct_RITs.rds')


