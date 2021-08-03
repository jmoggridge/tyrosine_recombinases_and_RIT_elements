
library(tidyverse)

rowAny <- function(x) rowSums(x) > 0 

# gets reverse complement of dna sequence
revcomp <- function(dna){
  map_chr(dna, ~{
    .x |> str_split('', simplify = T) |>
      seqinr::comp(ambiguous = T, forceToLower = F) |> 
      rev() |> 
      str_c(collapse = '') })
}

# Format and label seqs, reorient reverse complements, align
do_alignment <- function(df, seq, names, browse=FALSE, verbose = FALSE){
  seqs <- Biostrings::DNAStringSet(df |> pull({{seq}}))
  names(seqs) <- df |> pull({{names}})
  seqs <-  DECIPHER::OrientNucleotides(seqs, verbose = verbose)
  seqs.aligned <- DECIPHER::AlignSeqs(seqs, verbose = verbose)
  if (browse==TRUE){
    DECIPHER::BrowseSeqs(seqs.aligned)
  }
  return(seqs.aligned)
}


rit_elements <- read_rds('./results/rit_elements.rds')
glimpse(rit_elements)

rit_flanks <- 
  rit_elements |> 
  select(rit_id, matches('dna')) |> 
  rename_with(~str_remove(.x, 'rit_')) |> 
  filter(rowAny(across(.cols = contains('dna_'), ~!is.na(.x)))) |> 
  mutate(rc_downstream = map_chr(dna_downstream, revcomp)) |> 
  mutate(across(.cols = contains('dna_'), 
                .fns = ~map_dbl(.x, nchar), 
                .names = '{.col}_len')) |> 
  filter(dna_upstream_len > 100 & dna_downstream_len >100)

rit_flanks  |> 
  select(contains('len')) |> 
  summary()


# dna_aln <- do_alignment(rit_elements, rit_dna, rit_id, browse = T, verbose = T)
# 
# rit_elements |> 
#   slice(313:322) |> 
#   select(nuc_id)





rm(rit_elements)

downstream_aln <- do_alignment(rit_flanks, rit_dna_downstream, rit_id,
                               browse = T, verbose = T)
upstream_aln <- do_alignment(rit_flanks, rit_dna_upstream, rit_id,
                             browse = T, verbose = T)
