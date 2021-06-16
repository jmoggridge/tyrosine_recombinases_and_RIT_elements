


refseq <- read_rds('./data/non_integrase_seqs/refseq_non_integrases_raw.rds')
glimpse(refseq)
rm(queries, apikey, get_Efasta, get_ncbi_ids)


# count total number of seqs
refseq |> select(-name, -fasta) |> unnest(df) |> nrow()
# any duplicated records? -> none
refseq |> select(-name, -fasta) |> unnest(df) |>
  duplicated() |> any()
# any duplicated protein names? (title = acc + desc); no
refseq |> select(-name, -fasta) |> unnest(df) |>
  pull(title) |> unique() |> length()
# there are duplicated sequences though (285)
refseq |> select(-name, -fasta) |> unnest(df) |>
  pull(seq) |> unique() |> length()

# Keep 22250 unique sequences of 22535 in dataframe; keeping one occurrence where duplicate sequences are present
refseq <- refseq |> 
  select(-fasta) |> 
  unnest(df) |>
  group_by(seq) |> 
  sample_n(1) |> 
  # tidy up protein names for joining with hmmer results
  mutate(prot_acc = trimws(str_extract(title, '^.*? ')),
         prot_name = str_remove(title, '^.*? '))

# |> 
  # re-nest data and join to names,
  # nest(df = c(prot_acc, prot_name, title, seq)) 
# |> 
  # left_join(refseq |> select(name, fasta),
            # by = 'name')

glimpse(refseq)



#### Eukaryotic proteomes from Uniprot (non-integrases)

# Human, yeast and arabidopsis reference proteomes will be used as representative eukaryote sequences to classify as non-integrases. Reference proteomes have one sequence per protein.

# read in eukaryotic proteomes, then pfam non_integrases - seqs associated with other domains
set.seed(1)

other_non_integrases <- 
  tibble(
    path = Sys.glob('./data/non_integrase_seqs/uniprot_eukaryote_proteomes/*.fasta')) |> 
  mutate(name = str_remove_all(path, './.*/.*/|_proteome|.fasta')) |> 
  bind_rows(
    tibble(path = Sys.glob('./data/non_integrase_seqs/pfam_non_integrases/*.fa')) |> 
      mutate(name = str_remove_all(path, './.*/.*/|_[:alnum:]+_full|.fa'))
  ) |> 
  # extract titles and sequences
  mutate(fasta = map(path, ~readAAStringSet(.x))) |> 
  mutate(title = map(fasta, names),
         seq = map(fasta, paste)
  ) |> 
  select(-fasta) |> 
  unnest(c(title, seq)) |> 
  group_by(name) |> 
  sample_n(100) |> 
  ungroup()
  mutate(prot_name = paste0(name, '_', row_number()),
         prot_desc) |> 
  relocate(name, title, seq)


other_non_integrases |> 
  
glimpse(other_non_integrases) 
other_non_integrases |> 


non_integrases <- 
  read_rds('./data/non_integrases_df.rds')
non_integrases |> 
  ungroup() |> 
  
  slice(1:1000) |> 
  select(-fasta)

glimpse(non_integrases)

bind_rows(non_integrases, other_non_integrases)
# join refseq, uniprot, and pfam non_integrases

# remove duplicates

# downsample to 1000

# should remove this sequence: A0A1K1NN93_9FIRM, has a phage integrase domain




# combine stringsets and write a fasta file
non_integrases |> 
  pull(fasta) |> 
  reduce(c) |> 
  writeXStringSet('./data/all_non_integrase.fa')



# Write non-integrase df to file for later
write_rds(non_integrases, './data/non_integrases_df.rds')
rm(non_integrases, n)

