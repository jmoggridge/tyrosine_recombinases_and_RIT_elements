
library(seqinr)
library(tidyverse)
library(glue)

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

hamming_distance_lists <- function(x,y){
  if (length(x) != length(y)){
    message('Hamming distance: x and y not equal lengths')
    return(NA)
  }
  sum(x != y)
}
# hamming_distance('ABCD', 'ABBCD')
# hamming_distance('ABCD', 'BBCD') == 1

seed_locate <- function(seed, down){
  seed_size <- nchar(seed)
  sites <- map_lgl(
    .x = seq(nchar(down) - nchar(seed) + 1), 
    .f = ~str_sub(down, .x, .x + seed_size - 1) == seed
  )
  locations <- which(sites == T)
  if(!any(locations)) return(NA)
  else return(locations)
}

# filtering function for filter(rowAny(across(...)))
rowAny <- function(x) rowSums(x) > 0 




#### Data --------------------------------------------------------



## Active elements for terminal repeats ------------------------------------

nr_rits <- read_rds('results/non_redund_rits_NJ_clustered.rds')

# What is the copy number of distinct RIT coding sequences within a sequence record?
# Select nucleotides with > 1 copy; easiest and will have enough diversity of examples?
nr_rits |> 
  group_by(nuc_id, rit_id) |>  # identical sequences 
  mutate(rit_count = length(rit_dna)) |> 
  ungroup() |> 
  dplyr::count(rit_count) |> 
  arrange(desc(n))


# 1243 RIT dna sequences are distinct.
# There are 310 copies of 122 distinct (RIT dna, nucleotide record) pairs with >1 copies.
# Within this, there are 102 distinct RIT dna coding sequences.
# keep rits with multiple copies within sequence record.
multiple_rit_copies <- 
  nr_rits |> 
  # get first and last 100 bp of coding dna seq
  mutate(
    dna_first100 = map_chr(rit_dna, ~str_sub(.x, 1, 100)),
    dna_last100 = map2_chr(rit_dna, rit_length, ~str_sub(.x, .y-99, .y))
         ) |> 
  group_by(nuc_id, rit_id) |> 
  # identical coding sequences could have different upstream and downstream regions?
  mutate(rit_count = length(rit_dna)) |> 
  filter(rit_count > 1) |> 
  select(rit_count, rit_id, nuc_id, rit_dna_upstream, rit_dna_downstream,
         dna_first100, dna_last100) |> 
  distinct() |> 
  mutate(
    head = glue('{rit_dna_upstream}{dna_first100}'),
    tail = glue('{revcomp(rit_dna_downstream)}{revcomp(dna_last100)}'),
    head_tail = glue('{head}{tail}'),
    head_tail_name = glue('{rit_id}_{nuc_id}')
  ) |> 
  ungroup() |> 
  filter(nchar(head_tail) == 600) |> 
  arrange(rit_id, nuc_id)

## apply terminal inverted repeat finder to these active RITs
write_rds(multiple_rit_copies, './data/multiple_rit_copies.rds')


# Multiple copies of same RIT coding sequence within a single sequence record (contig, scaffold, or genome - depending on assembly).
# Removed any missing upstream or downstream sequence (4).
multiple_rit_copies <- 
  read_rds('./data/multiple_rit_copies.rds')
glimpse(multiple_rit_copies)


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
DECIPHER::BrowseSeqs(fasta)

# ## my data
# nr_rits <- read_rds('results/non_redundant_rits.rds')
# glimpse(nr_rits)
# 
# rit_flanks <- 
#   nr_rits |> 
#   transmute(rit_id, nuc_id, rit_dna, 
#             dna_up = rit_dna_upstream,
#             dna_down = rit_dna_downstream) |> 
#   distinct() |> 
#   filter(rowAny(across(.cols = contains('dna_'), ~!is.na(.x)))) |> 
#   mutate(rc_down = map_chr(dna_down, revcomp)) |> 
#   mutate(across(.cols = contains('dna_'), 
#                 .fns = ~map_dbl(.x, nchar), 
#                 .names = '{.col}_len')) |> 
#   filter(dna_up_len > 100 & dna_down_len >100)
# 
# rit_flanks  |> 
#   select(contains('len')) |> 
#   summary()
# 
# rm(rit_elements)






# 'imperfect palindromic repeats'

## NICOLE's data from paper

# Bolded bases are direct repeats contained within the terminal inverted repeats and for which there is an inverted copy at a precise distance in the direction of the recombinase genes.
strain_reps <- tribble(
  ~strain, ~five_prime_sequence,
  "Burkholderia phytofirmans Olga172",
  "TTATGCCGATTCCCGGATTATGCCG",
  "Cupriavidus metallidurans CH34", 
  "TTATGCCGACTCCCCGATTATGCCG",
  "Burkholderia sp. Ch1-1",
  "TTATGCCGACTTCCCGATTATGCCG",
  "Caulobacter sp. K31",
  "TAATGCCGCGATCCGGATTATGCCG",
  "Acidiphillium multivorum AIU301",
  "TAATGCCGAGATCCGGATTATGCCG",
  "Bifidobacterium longum NCC2705",
  "TTAAGCCGGGTTTGTTGTTAAGCCG",
  "Frankia sp. EANpec1",
  "TTATGCCGAGGGCCGGGTTATGCCG",
  "Novosphingobium sp. PP1Y",
  "TAATGCCGTGACCCGGATTATGCCG",
  "Candidatus S. usitatus Ellin6076",
  "ACTATGCCGCGTCCCGGACTATGCCGCGT",
  "Gramella forsetii KT0803 5’ Sequence",
  "ATTATGTAAAGTAAATTATTATGTAAAGT"
)

#### --------------------------------------------------------


# up_ans <- 'TTTGATTATGACAAGTAATATGTTATGTAAAGTGCAAG'
# down_ans <- 'TTCGATTATGACAAGTAATATGTTATGTAAAGTGCGAG'

# up <- rit_flanks$dna_up[[2]]
# down <- revcomp(rit_flanks$dna_down[[2]])

up <- 
  "TTTTTCTCTCCATGCCAACCCTCGAATTTATGAGTAGCTTATCATAATTTTGATTATGACAAGTAATATGTTATGTAAAGTGCAAGACTTATAAAGTGAGCAAACTTATTGAATATTAAAGAATTTATTTAGTCCCTTCGTATAAAAACTTCTCATTATATTGCATTCTGTATATATCCTAATTATTAAAGAATGATTTTA"

down <- 
  "CCGCTGCGGTCTTTTTTTCAGGTATGCCTTCAAGGGCGGAGAATCCGGTAAAGGCAAGTCCGGACAAGTTCAGTCTTGCGGTCGAGTGTGTCAAGCGATTCGATTATGACAAGTAATATGTTATGTAAAGTGCGAGTAAAAGTTGGTTTATAATCAAGTTGTTATGAACCAACTTGTCACTTACACTTTACATAATATTAT"

up; down

# take seed from x and locate in y
# find seed from in y seq

# matches <- tibble(
#   x_start = integer(),
#   y_start = integer(),
#   length = integer(),
#   hamming = integer()
# )
# matches

# seed_matches <- 
#   tibble(
#     x = up,
#     y = down,
#     pos_x = map_dbl(seq(nchar(up) - seed_size + 1), ~.x),
#     seed = map_chr(
#       .x = seq(nchar(up) - seed_size + 1), 
#       .f = ~ str_sub(up, .x, .x+seed_size)
#     )
#   ) |> 
#   mutate(pos_y = map_dbl(
#     .x = seed,
#     .f = ~seed_locate(.x, down)
#   )) |> 
#   filter(!is.na(pos_y)) |> 
#   mutate(
#     continue = ifelse(lag(pos_y) == (pos_y -1), T, F),
#     continue = ifelse(is.na(continue), F, continue)
#     ) |> 
#   filter(!continue) |> 
#   select(x, y, seed, pos_x, pos_y)
# seed_matches

get_seed_matches <- function(x, y, seed_size){
  seed_matches <- 
    tibble(
      x = x,
      y = y,
      pos_x = map_dbl(seq(nchar(up) - seed_size + 1), ~.x),
      seed = map_chr(
        .x = seq(nchar(up) - seed_size + 1), 
        .f = ~ str_sub(up, .x, .x+seed_size)
      )
    ) |> 
    mutate(pos_y = map_dbl(
      .x = seed,
      .f = ~seed_locate(.x, down)
    )) |> 
    filter(!is.na(pos_y)) |> 
    mutate(
      continue = ifelse(lag(pos_y) == (pos_y -1), T, F),
      continue = ifelse(is.na(continue), F, continue)
    ) |> 
    filter(!continue) |> 
    select(x, y, seed, pos_x, pos_y)
  if (nrow(seed_matches) == 0) return(NA)
  else return(seed_matches)
}



# TODO fix extension
extend_perfect_match <- function(up, down, pos_x, pos_y){
  up_ls <- up |> str_sub(pos_x, -1L) |> str_split('') |> unlist()
  down_ls <- down |> str_sub(pos_y, -1L) |> str_split('') |> unlist()
  
  perfect_match_length <- 
    tibble(match_length = seq(1, min(length(up_ls), length(down_ls)))) |> 
    mutate(hamming = map_dbl(
      .x = match_length,
      .f = ~hamming_distance_lists(up_ls[1:.x], down_ls[1:.x])
    )) |> 
    filter(hamming == 0) |> 
    filter(match_length == max(match_length)) |> 
    pull(match_length)
  return(perfect_match_length)
}

get_seed_matches(up, down, 10) |> 
  mutate(perf_match = pmap_dbl(
    .l = list(x, y, pos_x, pos_y),
    .f = function(x, y, pos_x, pos_y)
      extend_perfect_match(x, y, pos_x, pos_y)
  ))

perfect_matches <- function(x, y, seed_size){
  get_seed_matches(x, y, 10) 
  # |> 
  #   mutate(perf_match = pmap_dbl(
  #     .l = list(x, y, pos_x, pos_y),
  #     .f = function(x, y, pos_x, pos_y)
  #       extend_perfect_match(x, y, pos_x, pos_y)
  #   ))
}



multiple_rit_copies |> 
  mutate(inverted_rep_perf = map2(rit_dna_up, rc_down, ~perfect_matches(.x,.y, 10)))
  


pos_x <- seed_matches$pos_x
pos_y <- seed_matches$pos_y

up_ls <- up |> str_sub(pos_x, -1L) |> str_split('') |> unlist()
down_ls <- down |> str_sub(pos_y, -1L) |> str_split('') |> unlist()

perfect_match_length <- 
  tibble(match_length = seq(seed_size, min(length(up_ls), length(down_ls)))) |> 
  mutate(hamming = map_dbl(
    .x = match_length,
    .f = ~hamming_distance_lists(up_ls[1:.x], down_ls[1:.x])
    )) |> 
  filter(hamming == 0) |> 
  filter(match_length == max(match_length)) |> 
  pull(match_length)


ham <- 0
while(ham < 3){
  x <- str_sub(up, )
}














