# Clean up smart data from two fasta files for each subfamily of integrase, one 
# with the full protein sequences and the other just their integrase domains.

# Notes: 
# 1. Some domains are assigned to multiple classes. Thus, to ensure that no overlap between classes was present, I removed all the proteins having multiple subfamily labels for simplicity.
# 2. Not all domain sequences are unique; all protein sequences are unique
# 3. The information in the two fasta files doesn't match beyond the accession identifier (acc). Besides the domain names having the position of the domain within the protein, some of the description text following this doesn't match between the domain and the parent protein. Note this can cause joining to fail unless the keys are stipulated as acc

library(Biostrings)
library(tidyverse)
library(janitor)

## Domains data -------------------------------------------------

# Make dataframe of fastas for each subfamily of domains
# Split fasta names and seqs
domains <- 
  tibble(path = Sys.glob('./data/SMART/domain_fasta/*.fasta')) |>
  mutate(
    subfamily = str_remove_all(path, './data/SMART/domain_fasta/|\\.fasta'),
    dom_fasta = map(path, readAAStringSet),
    dom_n = map_int(dom_fasta, ~length(.x)),
    dom_name = map(dom_fasta, ~names(.x)),
    dom_seq = map(dom_fasta, paste)
  ) |> 
  select(-path, -dom_fasta) |> 
  select(subfamily, everything())
domains


# Unnest subfamilies and parse header information
domains <- domains |> 
  select(subfamily, dom_name, dom_seq) |>
  unnest(c(dom_name, dom_seq)) |> 
  mutate(dom_head = str_split(dom_name, '\\|'),
         acc = map_chr(dom_head, 3),
         acc_alt = map_chr(dom_head, 4) |> str_extract('^[^/]+'),
         desc = map_chr(dom_head, 4) |> str_remove('^.*?/'),
         start = str_extract(desc, '^[0-9]+'),
         stop = str_remove(desc, '^[0-9]+-') |> str_extract('[0-9]+'),
         description = str_remove(desc, '.*?[0-9]+-[0-9]+\\s')) |> 
  select(-c(dom_head, desc)) |> 
  relocate(subfamily, contains('acc'), start, stop, description, everything())

# domains |> View()

## # Check uniqueness: multiple class assignment per protein was found
## nrow(domains)
## length(unique(domains$dom_seq)) 
## length(unique(domains$acc))
## length(unique(domains$acc_alt))

# list of accessions with > 1 annotation to discard
garbage_annotations <- domains |> 
  group_by(acc) |> 
  filter(n() > 1) |> 
  ungroup()

# remove those proteins with multiple annotations
filtered_domains <- domains |> 
  anti_join(garbage_annotations)
filtered_domains

nrow(filtered_domains)
length(unique(filtered_domains$acc))
length(unique(filtered_domains$dom_seq))

# just to note:
# there are only 114,032 unique seqs (of 114,848) & some have multiple parent proteins
filtered_domains |> 
  group_by(dom_seq) |> 
  slice(1) |> 
  ungroup()

# some domain seqs belong to as many 8 different proteins
filtered_domains |>
  group_by(dom_seq) |>
  summarize(acc = length(unique(acc))) |>
  filter(acc >1) |> 
  arrange(desc(acc))

rm(domains)

## Protein data -------------------------------------------------

# Read the full protein sequences for 20 subfamilies from SMART
# Split fasta names and seqs
proteins <- 
  tibble(
    path = Sys.glob('./data/SMART/full_protein_fasta/*.fasta')) |> 
  mutate(
    subfamily = str_remove_all(
      path, './data/SMART/full_protein_fasta/|\\_proteins.fasta'),
    prot_fasta = map(path, readAAStringSet),
    prot_n = map_int(prot_fasta, ~length(.x)),
    prot_name = map(prot_fasta, ~names(.x)),
    prot_seq = map(prot_fasta, paste)
  ) |> 
  select(-path, -prot_fasta) |> 
  relocate(subfamily, everything())

proteins

# Parse header info from fastas
proteins <- proteins |> 
  select(-prot_n) |> 
  unnest(cols = c(prot_name, prot_seq)) |> 
  mutate(head = str_split(prot_name, '\\|'),
         acc = map_chr(head, 2),
         acc2 = map_chr(head, 3),
         acc_alt = str_extract(acc2, '^[^ ]+'),
         description = str_remove(acc2, '^[^ ]+') |> str_squish(),
         ) |> 
  select(-head, -acc2) |> 
  relocate(subfamily, contains('acc'), description) 

proteins

# remove proteins whose domains assigned to multiple subfamilies
filtered_proteins <- proteins |>
  anti_join(garbage_annotations,
            by = c("subfamily", "acc"))
filtered_proteins



## Join prot & domain ---------------------------------------------

# join filtered domains with their parent protein...
# but! some accessions don't match perfectly...
smart_df <- filtered_proteins |>
  left_join(filtered_domains,
            by = c("subfamily", "acc", "acc_alt","description"))

glimpse(smart_df)
write_rds(smart_df, './data/SMART/smart_df.rds', compress = 'gz')


## # All protein accessions and sequences are unique!
## smart_df |> pull(acc) |> unique() |> length()
## smart_df |> pull(prot_seq) |> unique() |> length()
## # Rows matched up correctly and no missing values created
## nrow(smart_df) == nrow(filtered_domains)
## nrow(smart_df) ==  nrow(filtered_proteins) 
## sum(is.na(smart_df))
## 

rm(proteins, filtered_domains, filtered_proteins, garbage_annotations)



## Split data ------------------------------------------------

# Might as well split data now

library(rsample)

set.seed(54321)
df_split <- initial_split(smart_df, 0.75, strata = subfamily)

smart_train <- training(df_split)
write_rds(smart_train, './data/SMART/smart_train.rds', compress = 'gz')

smart_test <- testing(df_split)
write_rds(smart_test, './data/SMART/smart_test.rds', compress = 'gz')

rm(df_split)


## Write fasta files for each domain subfamily ----------------
make_domain_fasta <- 
  function(df){
    fasta <- Biostrings::AAStringSet(df |> pull(dom_seq))
    names(fasta) <- df$acc
    dest <- paste0('./data/SMART/training_domain_fasta/', df$label[1], '.fa')
    writeXStringSet(fasta, filepath = dest)
  }

# create set fasta files of training domains for alignment
system('mkdir ./data/SMART/training_domain_fasta')
smart_train |> 
  mutate(label = subfamily) |> 
  group_by(subfamily) |> 
  group_walk(.f = ~make_domain_fasta(.x))
         

## train and test proteins fasta --------------------------------
make_protein_fasta <- 
  function(df, outfile){
    fasta <- Biostrings::AAStringSet(df |> pull(dom_seq))
    names(fasta) <- df$acc
    Biostrings::writeXStringSet(fasta, filepath = outfile)
  }

smart_train |> make_protein_fasta(outfile = './data/SMART/smart_train.fa')
smart_test |> make_protein_fasta(outfile = './data/SMART/smart_test.fa')


#### proceed to alignment script.... ----

