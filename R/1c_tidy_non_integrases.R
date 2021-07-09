### Script to combine non-integrase sequences from refseq, pfam, and uniprot, to provide negative examples for the classifier.

library(tidyverse)
library(Biostrings)

# Refseq dataset -----------------------------------------------------

refseq <- 
  # open up refseq raw data
  read_rds('./data/non_integrase_seqs/refseq_non_integrases_raw.rds') |> 
  mutate(group = name) |> 
  select(-name) |> 
  select(-fasta) |> 
  unnest(df) |>
  # remove duplicate sequences, keep first occurences
  group_by(seq) |> 
  slice_max(1) |> 
  ungroup() |> 
  # tidy up protein names for joining with hmmer results
  mutate(
    prot_name = str_extract(title, '^.*? ') |> str_squish(),
    prot_description = str_remove(title, '^.*? '),
    subfamily = 'other - refseq'
         )

glimpse(refseq)

## Pfam and Proteomes --------------------------------------------------

# Eukaryotic proteomes from Uniprot (non-integrases) and Pfam families 

# Human, yeast and arabidopsis reference proteomes will be used as representative eukaryote sequences to classify as non-integrases. Reference proteomes have one sequence per protein.

# read in fasta files: eukaryotic proteomes, then pfam non_integrase families

# pfam fastas
pfam <- 
  tibble(path = Sys.glob('./data/non_int*/pfam_non_int*/*.fa')) |> 
  mutate(
    group = str_remove_all(path, './.*/.*/|_full|.fa') |> str_squish(),
    subfamily = 'other - pfam: transposes, etc' 
    )

pfam

# uniprot proteome fastas
proteomes <- 
  tibble(path = Sys.glob('./data/non_int*/uniprot_*/*fasta')) |> 
  mutate(
    group = str_remove_all(path, './.*/.*/|_proteome|.fasta'),
    subfamily = 'other - eukaryote proteomes'
    )

proteomes

set.seed(1234)
# combine pfam and proteomes; read fastas; extract titles and sequences
other_non_integrases <- 
  bind_rows(pfam, proteomes) |> 
  mutate(fasta = map(path, ~readAAStringSet(.x))) |>
  mutate(title = map(fasta, names),
         seq = map(fasta, paste)
  ) |> 
  unnest(c(title, seq)) |> 
  select(-path, -fasta) |> 
  # should remove this sequence, it has a phage integrase domain
  filter(!str_detect(title, 'A0A1K1NN93_9FIRM')) |> 
  # downsample to 2000 per group
  group_by(group) |>
  sample_n(size = 2000) |>
  ungroup() |>
  # create accession keys as '<family>_<row#>'
  mutate(prot_name = paste0(group, '_', row_number()),
         prot_description = title)

other_non_integrases
glimpse(other_non_integrases)
rm(pfam, proteomes)

## Join ------------------------------------------------------------

# join all the data; match names with integrases dataframe
non_integrases <- 
  bind_rows(refseq, other_non_integrases) |> 
  transmute(subfamily,
            group, 
            acc = prot_name, 
            prot_seq = seq |> str_to_upper(), 
            title, 
            description = prot_description) |> 
  # remove any with ambiguous bases
  filter(!str_detect(prot_seq, 'B|J|O|U|X|Z')) |> 
  # remove any that are very short
  filter(!nchar(prot_seq) < 30) |> 
  # keep only unique sequences; 
  group_by(prot_seq) |> 
  sample_n(1) |> 
  ungroup()

glimpse(non_integrases)
rm(refseq, other_non_integrases)


## Data checks ------------------------------------------------------------

# verify id and sequences are unique
non_integrases |> distinct() |> nrow()
non_integrases |> pull(acc) |> unique() |> length()
non_integrases |> pull(prot_seq) |> unique() |> length()

# check group sizes
non_integrases |> 
  count(group, sort = T) |> 
  print.data.frame()

# check subfamily sizes
non_integrases |> 
  count(subfamily, sort = T) |> 
  print.data.frame()

# plot seq lengths
non_integrases |> 
  mutate(length = nchar(prot_seq)) |> 
  ggplot(aes(length, fct_rev(group))) +
  geom_jitter(alpha = 0.2, size = 0.2) +
  scale_x_log10() +
  theme_bw() +
  labs(y = '', x = 'length (residues)', subtitle = 'Non-integrase protein lengths')

non_integrases |> 
  mutate(length = nchar(prot_seq)) |> 
  ggplot(aes(length, fct_rev(subfamily))) +
  geom_jitter(alpha = 0.2, size = 0.2) +
  scale_x_log10() +
  theme_bw() +
  labs(y = '', x = 'length (residues)', subtitle = 'Non-integrase protein lengths')

# what sequence is 6k residues long?
non_integrases |> 
  mutate(length = nchar(prot_seq)) |> 
  filter(length == max(length)) |> 
  pull(title)

# verify no duplicate sequences
length(unique(non_integrases$prot_seq)) == nrow(non_integrases)

write_rds(non_integrases, './data/non_integrase_seqs/nonint_df.rds')

