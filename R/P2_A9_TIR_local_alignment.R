
library(tidyverse)
library(glue)
library(Biostrings)
library(furrr)

# gets reverse complement of dna sequence
revcomp <- function(dna){
  map_chr(dna, ~{
    .x |> str_split('', simplify = T) |>
      seqinr::comp(ambiguous = T, forceToLower = F) |> 
      rev() |> 
      str_c(collapse = '') })
}


## Active elements for terminal repeats ------------------------------------

nr_rits <- read_rds('results/non_redund_rits_NJ_clustered.rds') |> 
  mutate(obs_id = row_number())

# What is the copy number of distinct RIT coding sequences within a sequence record?
# Select nucleotides with > 1 exact copy of a rit coding sequence;
# should be easiest and will have enough diversity of examples for inverted repeats?
nr_rits |> 
  group_by(nuc_id, rit_id) |>  # identical sequences 
  mutate(rit_count = length(rit_dna)) |> 
  ungroup() |> 
  dplyr::count(rit_count) |> 
  arrange(desc(n))
nr_rits

# 915 RIT dna sequences are distinct.
nr_rits |> dplyr::count(rit_dna) |> nrow()

# 1364 distinct rit_dna, nuc_id pairs
nr_rits |> dplyr::count(rit_dna, nuc_id) |>  nrow()
# There are 122 distinct (RIT dna, nucleotide record) pairs with >1 copies.
nr_rits |> dplyr::count(rit_dna, nuc_id) |> filter(n>1) |>  nrow()
# There are 317 copies total in this cohort
nr_rits |> dplyr::count(rit_dna, nuc_id) |> filter(n>1) |>  pull(n) |> sum()
# Of 104 distinct RIT dna coding sequences.
nr_rits |> dplyr::count(rit_dna, nuc_id) |> filter(n>1) |> dplyr::count(rit_dna, name = 'nn')




## apply terminal inverted repeat finder to these active RITs --------

# keep rits with multiple copies within sequence record
# Keep the sequence records with multiple copies of a RIT coding sequence within a single  (contig, scaffold, or genome - depending on assembly). Exact matches only (recently mobile).
# Removed any missing upstream or downstream sequence (5 obs).
multiple_rit_copies <- 
  nr_rits |> 
  # first and last 100 bp of coding dna seq
  mutate(
    dna_first100 = map_chr(rit_dna, ~str_sub(.x, 1, 100)),
    dna_last100 = map2_chr(rit_dna, rit_length, ~str_sub(.x, .y-99, .y))
  ) |> 
  group_by(nuc_id, rit_id) |> 
  # identical coding sequences could have different upstream and downstream regions?
  mutate(rit_count = length(rit_dna)) |> 
  filter(rit_count > 1) |> 
  select(rit_count, rit_id, nuc_id, obs_id,
         rit_dna_upstream, rit_dna_downstream,
         dna_first100, dna_last100) |> 
  distinct() |> 
  mutate(
    head = glue('{rit_dna_upstream}{dna_first100}'),
    tail = glue('{revcomp(rit_dna_downstream)}{revcomp(dna_last100)}'),
    head_tail = glue('{head}{tail}'),
    head_tail_name = glue('{rit_id}_{nuc_id}')
  ) |> 
  ungroup() |> 
  filter(nchar(head_tail) == 1200) |> 
  arrange(nuc_id, rit_id, obs_id)



## ** The active elements aren't all the same clusters! ** -----
## 102 active elements from 94 clusters @ dist 0.005; 
multiple_rit_copies |> dplyr::count(rit_id) |> nrow()

multiple_rit_copies |> 
  left_join(nr_rits |> select(rit_id, cluster0_005NJ) |> distinct()) |> 
  dplyr::count(cluster0_005NJ) |> 
  arrange(desc(n))

# 93 clusters @ 0.01
multiple_rit_copies |> 
  left_join(nr_rits |> select(rit_id, cluster0_01NJ) |> distinct()) |> 
  dplyr::count(cluster0_01NJ) |> 
  arrange(desc(n))
# 72 @ 0.3
multiple_rit_copies |> 
  left_join(nr_rits |> select(rit_id, cluster0_3NJ) |> distinct()) |> 
  dplyr::count(cluster0_3NJ) |> 
  arrange(desc(n))

# multiple_rit_copies |> View()
skimr::skim(multiple_rit_copies)

# create forward and rev complement seq to search
# provide -500: +100 from p1 start and -100:+500 from p3 end of RIT CDS
multiple_rit_copies <- 
  multiple_rit_copies |> 
  mutate(
    upstream = glue('{rit_dna_upstream}{dna_first100}'),
    rc_downstream = revcomp(glue('{dna_last100}{rit_dna_downstream}')),
  ) 

write_rds(multiple_rit_copies, './data/multiple_rit_copies.rds')



## examine sequences - can see where repeats stop and sequence isn't within some larger
# element
fasta <- Biostrings::DNAStringSet(multiple_rit_copies$head_tail)
names(fasta) <- multiple_rit_copies$head_tail_name
# DECIPHER::BrowseSeqs(fasta)


do_local_alignment <- function(x, y){
  x <- DNAString(x)
  y <- DNAString(y)
  aligned <-pairwiseAlignment(x, y, type = 'local', gapOpening = 4, gapExtension = 1)
  return(
    tibble(
      al_score = score(ali),
      al_percent_id = pid(aligned),
      al_length = nchar(aligned),
      al_match = nmatch(aligned),
      al_mismatch = nmismatch(aligned),
      al_edit = nedit(aligned),
      al_fwd_seq = paste(pattern(aligned)),
      al_rc_seq = paste(subject(aligned)),
      al_fwd_start = str_locate(x, toString(pattern(aligned)))[1,1],
      al_fwd_end = str_locate(x, toString(pattern(aligned)))[1,2],
      al_rc_start = str_locate(x, toString(subject(aligned)))[1,1],
      al_rc_end = str_locate(x, toString(subject(aligned)))[1,2],
      al_align = list(aligned)
    )
    )
}


rit_IRs <- multiple_rit_copies |> 
  select(obs_id, upstream, rc_downstream) |> 
  mutate(alignment = future_map2(upstream, rc_downstream, do_local_alignment)) |> 
  unnest_wider(col = alignment) |> 
  right_join(multiple_rit_copies)

glimpse(rit_IRs)

rit_IRs_df <- nr_rits |> 
  right_join(rit_IRs) |> 
  mutate(
    ir_5p_start =  -500 + al_fwd_start,
    ir_5p_end = -500 + al_fwd_end,
    ir_3p_start =  500 - al_rc_start,
    ir_3p_end = 500 - al_rc_end,
  )

write_rds(rit_IRs_df, 'results/rit_IRs.rds')

# length, matches
rit_IRs_df |> 
  pivot_longer(c(al_match, al_percent_id, al_length)) |> 
  ggplot(aes(value)) +
  geom_histogram() +
  facet_wrap(~name, scales = 'free') +
  theme_light() +
  labs(title = 'Inverted repeat alignments for active RITs (multiple-copy per contig)')

rit_IRs_df |> 
  pivot_longer(ir_5p_start:ir_3p_end) |> 
  mutate(
    side = ifelse(str_detect(name, '5p'), 'upstream', 'downstream'),
    startstop = ifelse(str_detect(name, 'start'), 'start', 'stop')) |> 
  
  ggplot(aes(value, fill = startstop, color = startstop)) +
  geom_density(alpha =0.2) +
  facet_wrap(~fct_rev(side), scales = 'free') +
  theme_light() +
  labs(
    subtitle =
      'Where are the aligned inverted repeats located relative to beginning and end \
of the coding sequence?', 
    color = NULL, fill = NULL)


# show some sample alignments....

rit_IRs_df$al_align[[1]]



# 
# x <- DNAString(df$upstream[101])
# y <- DNAString(df$rc_downstream[101])
# 
# ali <- pairwiseAlignment(x, y, type = 'local', gapOpening = 3, gapExtension = 2)
# 
# Views(ali)
# aligned(ali)
# score(ali)
# summary(ali)
# consensusString(ali, )
# 
