
library(tidyverse)

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

# Active elements!
# multiple copies and not missing upstream or downstream
multiple_rit_copies <- 
  read_rds('./data/multiple_rit_copies.rds') |> 
  filter(nchar(head_tail) == 400) |> 
  ungroup() |> 
  arrange(rit_id)

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














