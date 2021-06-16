### Script to combine non-integrase sequences from refseq, pfam, and uniprot, to provide negative examples for the classifier.

library(tidyverse)
library(Biostrings)


refseq <- 
  # open up refseq raw data
  read_rds('./data/non_integrase_seqs/refseq_non_integrases_raw.rds') |> 
  select(-fasta) |> 
  unnest(df) |>
  # remove duplicate sequences, keep first occurences
  group_by(seq) |> 
  slice_max(1) |> 
  ungroup() |> 
  # tidy up protein names for joining with hmmer results
  mutate(prot_acc = trimws(str_extract(title, '^.*? ')),
         prot_name = str_remove(title, '^.*? '))

glimpse(refseq)


#### Eukaryotic proteomes from Uniprot (non-integrases) and Pfam families 

# Human, yeast and arabidopsis reference proteomes will be used as representative eukaryote sequences to classify as non-integrases. Reference proteomes have one sequence per protein.

# read in fasta files: eukaryotic proteomes, then pfam non_integrase families
set.seed(1)

# pfam fastas
pfam <- 
  tibble(path = Sys.glob('./data/non_int*/pfam_non_int*/*.fa')) |> 
  mutate(name = str_remove_all(path, './.*/.*/|_full|.fa'))

# uniprot proteome fastas
proteomes <- 
  tibble(path = Sys.glob('./data/non_int*/uniprot_*/*fasta')) |> 
  mutate(name = str_remove_all(path, './.*/.*/|_proteome|.fasta'))

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
  # downsample to 1000 per group
  group_by(name) |> 
  sample_n(2000) |> 
  ungroup() |> 
  # create accession keys as '<family>_<row#>'
  mutate(prot_acc = paste0(name, '_', row_number()),
         prot_name = title)

# join all the data
non_integrases <- 
  bind_rows(refseq, other_non_integrases) |> 
  group_by(seq) |> 
  sample_n(1) |> 
  ungroup() |> 
  relocate(name, prot_acc, prot_name)

non_integrases |> 
  mutate(length = nchar(seq)) |> 
  ggplot(aes(length, fct_rev(name))) +
  geom_jitter(alpha = 0.2, size = 0.2) +
  labs(y = '', x = 'length (residues)', subtitle = 'Non-integrase protein sizes')

non_integrases |> 
  mutate(length = nchar(seq)) |> 
  filter(length == max(length)) |> 
  pull(prot_acc)

# verify no duplicate sequences
length(unique(non_integrases$seq)) == nrow(non_integrases)


# Write non-integrase df to file for later
write_rds(non_integrases, './data/non_integrases_df.rds')


non_int_fa <- AAStringSet(toupper(non_integrases$seq))
names(non_int_fa) <- non_integrases$prot_acc
non_int_fa
writeXStringSet(non_int_fa, './data/non_integrases.fa')

rm(refseq, pfam, proteomes, other_non_integrases)



