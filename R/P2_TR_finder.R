
### Setup ------------------------------------------------------------

library(seqinr)
library(tidyverse)
library(glue)
library(beepr)
library(progress)

#### Functions --------------------------------------------------------

# gets reverse complement of dna sequence
revcomp <- function(dna){
  map_chr(dna, ~{
    .x |> str_split('', simplify = T) |>
      seqinr::comp(ambiguous = T, forceToLower = F) |> 
      rev() |> 
      str_c(collapse = '') })
}

# number of mismatches in two equal length strings
hamming_distance <- function(x,y){
  if (nchar(x) != nchar(y)){
    message('Hamming distance: x and y not equal lengths')
    return(NA)
  }
  sum(unlist(str_split(x, '')) != unlist(str_split(y, '')))
}

# hamming distance for 2 eq. length lists of characters
hamming_distance_lists <- function(x,y){
  if (length(x) != length(y)){
    message('Hamming distance: x and y not equal lengths')
    return(NA)
  }
  sum(x != y)
}
# hamming_distance('ABCD', 'ABBCD')
# hamming_distance('ABCD', 'BBCD') == 1


# for given seed from fwd seq, locate exact match position in rev seq.
seed_locate <- function(seed, revr){
  seed_size <- nchar(seed)
  sites <- map_lgl(
    .x = seq(nchar(revr) - nchar(seed) + 1), 
    .f = ~str_sub(revr, .x, .x + seed_size - 1) == seed
  )
  locations <- which(sites == T)
  if(!any(locations)) return(NA)
  else return(locations)
}

# filtering function for filter(rowAny(across(...)))
rowAny <- function(x) rowSums(x) > 0 



# fwd is sense strand, rev is antisense relative to RIT.
# ie. seq of interest and the rev complement seq to match
# return longest perfect IR
get_seed_matches <- function(fwd, rev, seed_size){

  seed_matches <- 
    # create list of position: kmer
    tibble(
      f = fwd,
      r = rev,
      pos_f = map_dbl(seq(nchar(fwd) - seed_size + 1), ~.x),
      seed = map_chr(
        .x = seq(nchar(fwd) - seed_size + 1), 
        .f = ~ str_sub(fwd, .x, .x + seed_size - 1)
      )
    ) |> 
    # find positions in y with exact match for each seed
    mutate(
      pos_r = map(
        .x = seed,
        .f = ~seed_locate(.x, rev)
      )) |> 
    unnest(pos_r) |> 
    # remove any seeds that failed to match
    filter(!is.na(pos_r)) |> 
    
    # continue == T is an extension or the previous row; 
    # in both x and y
    mutate(
      continue = ifelse(lag(pos_r) == (pos_r -1), T, F),
      continue = ifelse(is.na(continue), F, continue)
    ) |> 
    # only keep first position of the match
    filter(!continue) |> 
    select(f, r, seed, pos_f, pos_r)
  
  # if no match
  if (nrow(seed_matches) == 0)
    return(NA)
  
  # extend all matches; keep longest perfect match btwn forward & reverse
  # 
  else 
    best_match <- 
      seed_matches |> 
      mutate(match_length = pmap_dbl(
        .l = list(f, r, pos_f, pos_r),
        .f = function(f, r, pos_f, pos_r)
          extend_perfect_match(f, r, pos_f, pos_r)
      )) |> 
      filter(match_length == max(match_length)) |> 
      select(pos_f, pos_r, match_length) |> 
      mutate(ir_seq = str_sub(fwd, pos_f, pos_f + match_length - 1))
  
  return(best_match)
}




#
extend_perfect_match <- function(fwd, rev, pos_x, pos_y){
  
  x_ls <- fwd |> str_sub(pos_x, -1L) |> str_split('') |> unlist()
  y_ls <- rev |> str_sub(pos_y, -1L) |> str_split('') |> unlist()
  
  # get the perfect longest match fo
  perfect_match_length <-
    tibble(match_length = seq(1, min(length(x_ls), length(y_ls)))) |>
    mutate(hamming = map_dbl(
      .x = match_length,
      .f = ~hamming_distance_lists(x_ls[1:.x], y_ls[1:.x])
    )) |>
    filter(hamming == 0) |>
    filter(match_length == max(match_length)) |>
    pull(match_length)
  return(perfect_match_length)
}








#### Data --------------------------------------------------------



## Active elements for terminal repeats ------------------------------------

nr_rits <- read_rds('results/non_redund_rits_NJ_clustered.rds')

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

# 1243 RIT dna sequences are distinct.
# There are 310 copies of 122 distinct (RIT dna, nucleotide record) pairs with >1 copies.
# Within this, there are 102 distinct RIT dna coding sequences.


## apply terminal inverted repeat finder to these active RITs --------

# keep rits with multiple copies within sequence record
# Keep the sequence records with multiple copies of a RIT coding sequence within a single  (contig, scaffold, or genome - depending on assembly). Exact matches only (recently mobile).
# Removed any missing upstream or downstream sequence (4).
multiple_rit_copies <- 
  nr_rits |> 
  # first and last 100 bp of coding dna seq
  mutate(
    obs_id = row_number(),
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
skimr::skim(mulitple_rit_copies)

write_rds(multiple_rit_copies, './data/multiple_rit_copies.rds')


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


## examine sequences - can see where repeats stop and sequence isn't within some larger
# element
fasta <- Biostrings::DNAStringSet(multiple_rit_copies$head_tail)
names(fasta) <- multiple_rit_copies$head_tail_name
# DECIPHER::BrowseSeqs(fasta)


# 
# # 'imperfect palindromic repeats'
# 
# ## NICOLE's data from paper
# 
# # Bolded bases are direct repeats contained within the terminal inverted repeats and for which there is an inverted copy at a precise distance in the direction of the recombinase genes.
# strain_reps <- tribble(
#   ~strain, ~five_prime_sequence,
#   "Burkholderia phytofirmans Olga172",
#   "TTATGCCGATTCCCGGATTATGCCG",
#   "Cupriavidus metallidurans CH34", 
#   "TTATGCCGACTCCCCGATTATGCCG",
#   "Burkholderia sp. Ch1-1",
#   "TTATGCCGACTTCCCGATTATGCCG",
#   "Caulobacter sp. K31",
#   "TAATGCCGCGATCCGGATTATGCCG",
#   "Acidiphillium multivorum AIU301",
#   "TAATGCCGAGATCCGGATTATGCCG",
#   "Bifidobacterium longum NCC2705",
#   "TTAAGCCGGGTTTGTTGTTAAGCCG",
#   "Frankia sp. EANpec1",
#   "TTATGCCGAGGGCCGGGTTATGCCG",
#   "Novosphingobium sp. PP1Y",
#   "TAATGCCGTGACCCGGATTATGCCG",
#   "Candidatus S. usitatus Ellin6076",
#   "ACTATGCCGCGTCCCGGACTATGCCGCGT",
#   "Gramella forsetii KT0803 5’ Sequence",
#   "ATTATGTAAAGTAAATTATTATGTAAAGT"
# )

#### --------------------------------------------------------

# do one example ---

# up_ans <- 'TTTGATTATGACAAGTAATATGTTATGTAAAGTGCAAG'
# down_ans <- 'TTCGATTATGACAAGTAATATGTTATGTAAAGTGCGAG'

# up <- rit_flanks$dna_up[[2]]
# down <- revcomp(rit_flanks$dna_down[[2]])

fwd <-
  "TTTTTCTCTCCATGCCAACCCTCGAATTTATGAGTAGCTTATCATAATTTTGATTATGACAAGTAATATGTTATGTAAAGTGCAAGACTTATAAAGTGAGCAAACTTATTGAATATTAAAGAATTTATTTAGTCCCTTCGTATAAAAACTTCTCATTATATTGCATTCTGTATATATCCTAATTATTAAAGAATGATTTTA"

rev <-
  "CCGCTGCGGTCTTTTTTTCAGGTATGCCTTCAAGGGCGGAGAATCCGGTAAAGGCAAGTCCGGACAAGTTCAGTCTTGCGGTCGAGTGTGTCAAGCGATTCGATTATGACAAGTAATATGTTATGTAAAGTGCGAGTAAAAGTTGGTTTATAATCAAGTTGTTATGAACCAACTTGTCACTTACACTTTACATAATATTAT"

fwd; rev











pb <- progress::progress_bar$new(
  format = "Finding repeats: [:bar] :percent eta: :eta",
  total = nrow(multiple_rit_copies),
  clear = FALSE,
)

# create forward and rev complement seq to search
# provide -500: +100 from p1 start and -100:+500 from p3 end of RIT CDS
multiple_rit_copies <- 
  multiple_rit_copies |> 
  mutate(
    upstream = glue('{rit_dna_upstream}{dna_first100}'),
    rc_downstream = revcomp(glue('{dna_last100}{rit_dna_downstream}')),
  ) 

# find longest perfect matches in upstream and rc of downstream for each RIT
rit_inverted_reps <- 
  multiple_rit_copies |> 
  select(obs_id, upstream, rc_downstream) |> 
  mutate(
    longest_perfect_IR = map2(
      .x = upstream, 
      .y = rc_downstream,
      .f = ~{
        pb$tick()
        return(get_seed_matches(fwd = .x, rev = .y, seed_size = 8))
      }
    )) |> 
  unnest(cols = c(longest_perfect_IR))

beepr::beep()

# perfect repeat length
rit_inverted_reps |> 
  ggplot(aes(match_length)) +
  geom_histogram() +
  scale_x_log10() +
   labs(title = 'Perfect inverted repeat match length')

# distance from coding start for IR start and stop positions
rit_IR_relative_pos <- rit_inverted_reps |> 
  left_join(multiple_rit_copies |> select(rit_id, nuc_id, obs_id)) |> 
  mutate(
    ir5_start =  -500 + pos_f,
    ir5_end = -500 + pos_f + match_length - 1,
    ir3_start =  500 - pos_r,
    ir3_end = 500 - (pos_r + match_length - 1),
  )

rit_IR_relative_pos |> 
  pivot_longer(ir5_start:ir3_end) |> 
  mutate(
    side = ifelse(str_detect(name, 'ir5'), 'upstream', 'downstream'),
    startstop = ifelse(str_detect(name, 'start'), 'start', 'stop')) |> 
  
  ggplot(aes(value, fill = startstop, color = startstop)) +
  geom_density(alpha =0.2) +
  facet_wrap(~fct_rev(side), scales = 'free') +
  labs(subtitle ='Where are the perfect inverted repeats located relative\n to beginning and end of coding sequence?', color = NULL, fill = NULL)





rit_inverted_reps |> 
  print(n=300)



fwd <- multiple_rit_copies$upstream[[1]]
rc <- multiple_rit_copies$rc_downstream[[1]]
fwd; rc

perf_ir <- get_seed_matches(fwd, rc, seed_size = 8)
perf_ir

perf_length <- perf_ir$match_length[[1]]
perf_length

f_start <- perf_ir$pos_f[[1]]
f_start

f_end <- perf_ir$pos_f[[1]] + perf_length -1
f_end

r_start <- perf_ir$pos_r[[1]]
r_start

r_end <- perf_ir$pos_r[[1]] + perf_length -1
r_end

r_end - r_start + 1

str_sub(fwd, f_start, f_end)
str_sub(rev, r_start, r_end)


fwd_ls <- str_split(fwd, '') |> as_vector()
fwd_ls
rc_ls <- str_split(rc, '') |> as_vector()

# backwards from 5' of match

left_fwd <- rev(fwd_ls[1:f_start - 1])
left_rc <- rev(rc_ls[1:r_start - 1]) 

right_fwd <- fwd_ls[f_end + 1: length(fwd_ls)]
right_rc <- rc_ls[r_end + 1: length(rc_ls)]

left_fwd; left_rc

ham <- 1
pos_x <- 2
pos_y <- 2







# try_left <- function(left_fwd, left_rc){
#   # skip over mismatch at first position
#   ham <- 1
#   while (ham < 2){
#     
#   }
# }

paste(left_fwd, collapse = '')
paste(left_rc, collapse = '')



# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
