## 1a_tidy_smart_data.R

## Briefing -----

# Cleans up smart data from two fasta files for each subfamily of integrase, one 
# with the full protein sequences and the other just their integrase domains.

# Notes: 

# 1. Some domains (5704) come from objects (acc ids 2844) are assigned to multiple subfamilies in SMART. Thus, to ensure that no overlap between classes was present, I removed all the proteins having multiple subfamily labels for simplicity. There are 114,848  observations with unique ids and single subfamily labels after filtering.

# 2. Not all domain sequences are unique; there are 114,032 unique sequences from the total of 114848. However, all protein sequences are unique but some domain seqs are found in as many 8 different proteins.

# 3. The information in the two fasta files doesn't match beyond the accession identifier (acc). Besides the domain names having the position of the domain within the protein, some of the description text following this doesn't match between the domain and the parent protein. Note this can cause joining to fail unless the keys are stipulated as acc


library(Biostrings)
library(tidyverse)
library(janitor)

## Domains data -------------------------------------------------

# Make dataframe of fastas for each subfamily of domains
# Split fasta names and seqs into separate columns
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


domains |> pull(dom_n) |> sum()

domains |> pull(dom_n) |> median()


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
nrow(domains)                     # 120552
length(unique(domains$dom_seq))   # 119694, some shared domain seqs
length(unique(domains$acc))       # 117692, even fewer unique ids
length(unique(domains$acc_alt))   # same # of ids

# nrow(domains)  - 

# list of accessions with > 1 annotation to discard
garbage_annotations <- domains |> 
  group_by(acc) |> 
  filter(n() > 1) |> 
  ungroup()
garbage_annotations |> nrow()                 # 5704 entries from ...
garbage_annotations |> count(acc) |> nrow()   # 2844 ids w multiple labels

# remove those proteins with multiple annotations; 114,848 obs w unique ids
filtered_domains <- domains |> 
  anti_join(garbage_annotations)

filtered_domains
filtered_domains$acc |> unique() |> length()
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
smart_df |> pull(acc) |> unique() |> length()
smart_df |> pull(prot_seq) |> unique() |> length()
# Rows matched up correctly and no missing values created
nrow(smart_df) == nrow(filtered_domains)
nrow(smart_df) ==  nrow(filtered_proteins)
sum(is.na(smart_df))


## End. Proceed to alignment script.... ----


ggplot(smart_df, aes(fct_rev(subfamily))) +
  geom_bar() +
  theme_classic() +
  coord_flip() +
  labs(x = NULL)


