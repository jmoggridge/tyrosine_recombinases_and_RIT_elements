## P2_A2_repair_nucleotide_data -----


library(tidyverse)
source('./R/P2_entrez_functions.R')

# ids dataset from P2_A1
id_data <- 
  read_rds('./data/CDD/cdd_prot_nuc_tax_ids.rds') |> 
  filter(!is.na(nuc_id))
id_data


# index of nucleotide ids for sequence files
nuc_ids <- unique(id_data$nuc_id)
n <- length(nuc_ids)
nuc_sets <- 
  tibble(id_set = map(seq(1, n, 1000), ~nuc_ids[.x: min(.x + 999, n)])) |> 
  mutate(set_number = row_number())
nuc_sets

summaries <- 
  read_rds('./data/CDD/nuc_summary.rds') |> 
  mutate(nuc_id = name) |> 
  select(-c(name, gi, uid, term)) |> # all same as nuc_id
  mutate(across(everything(), ~na_if(.x, '')))
glimpse(summaries)
# skimr::skim(summaries)

## Issues? Some EDA ----

df <- id_data |>
  select(nuc_id) |> 
  distinct() |> 
  left_join(summaries)

# some sequences removed, 2 superceded, 2 updated
df |> 
  count(comment)

# 205 'suppressed'; 2 replaced
df |> 
  count(status)

# updated: 1034665437, 1632804275
df |> 
  filter(!is.na(replacedby)) |> 
  count(nuc_id, replacedby, comment)

# most circular are probably complete genomes or maybe plasmids
df |> 
  count(topology)
#  topology     n
# 1 circular  2060
# 2 linear   13425
# 3 not-set     10

id_data |> 
  left_join(summaries |> select(nuc_id, genome)) |> 
  count(genome)
# 1 ""            7217
# 2 "chromosome"  1403
# 3 "genomic"     6147
# 4 "plasmid"      678
# 5  NA             50

# two seqs are RNA... should exclude. possibly duplicate nuc records for same protein with RNA and DNA?
id_data |> 
  left_join(summaries |> select(nuc_id, moltype)) |> 
  count(moltype)
id_data |> 
  left_join(summaries |> select(nuc_id, biomol)) |> 
  count(biomol)

# 
id_data |> 
left_join(summaries |> select(nuc_id, tech)) |> 
  count(tech)
id_data |> 
  left_join(summaries |> select(nuc_id, sourcedb)) |> 
  count(sourcedb)

# seq lengths seem okay
# id_data |> 
#   left_join(summaries |> transmute(nuc_id, kbp = as.numeric(slen)/1000)) |> 
#   select(nuc_id, kbp) |> 
#   distinct() |>
#   ggplot(aes(kbp)) +
#   geom_histogram() +
#   scale_x_log10() +
#   labs(title = glue('Lengths of distinct nt seq'))



## Fix Superceded Records------
# 1887514880
# 1867089197


### get replacement data from entrez ----
# replace two superceded sequences "1887514880" "1867089197"
# these don't have any new id given in replaced by
superceded <- 
  id_data |> 
  left_join(
    summaries |>
      select(nuc_id, caption, comment, replacedby)
  ) |> 
  filter(str_detect(comment, 'CP053749|CP063356')) |> 
  mutate(new_nuc_id = str_extract(comment, 'CP053749|CP063356')) |> 
  mutate(
    token = map(
      new_nuc_id, 
      ~entrez_search(db = 'nuccore', term = .x, use_history = T)
    ),
    fasta = map(
      token,
      ~entrez_fetch(db = 'nuccore', web_history = .x$web_history, 
                    rettype = 'fasta')
      ),
    nuc_summary = map(
      token,
      ~entrez_summary(db = 'nuccore', web_history = .x$web_history)
    )
    ) |> 
  mutate(
    nuc_ss = map(fasta, fasta_to_DNAss),
    nuc_name = map(nuc_ss, names),
    nuc_seq = map(nuc_ss, paste)
  ) |> 
  select(nuc_id, new_nuc_id, nuc_name, nuc_seq, nuc_summary)

print(superceded)

### fixed seqs: superceded_seq ----

# prepare superceded replacement sequences
superceded_seq <- superceded |> 
  select(nuc_id, new_nuc_id, nuc_name, nuc_seq)

# prepare superceded replacement summaries
superceded_sum <-
  superceded |> 
  select(nuc_id, nuc_summary) |> 
  unnest_wider(nuc_summary) |> 
  mutate(across(everything(), as.character)) |> 
  mutate(nuc_id = uid) |> 
  select(-c(gi, uid, term))
  
# replace superceded records in summaries.
summaries <- 
  summaries |> 
  filter(!nuc_id %in% c("1887514880","1867089197")) |> 
  bind_rows(superceded_sum)

### fix id_data ----
# prepare ids and get new tax ids for replacement records
superced_ids <- 
  superceded |> 
  select(nuc_id, nuc_summary) |> 
  left_join(id_data) |> 
  unnest_wider(nuc_summary) |> 
  transmute(nuc_id = uid, cdd_id, prot_id) |> 
  link_nuccore_taxonomy(id = nuc_id)

id_data <- 
  id_data |> 
  filter(!nuc_id %in% c("1887514880","1867089197")) |> 
  bind_rows(superced_ids)

### fix taxonomy ----
superced_tax <- 
  superced_ids |> 
  fetch_taxonomy(id = tax_id)
superced_tax

read_rds('./data/CDD/tax_data.rds') |> 
  bind_rows(superced_tax) |> 
  distinct() |> 
  write_rds('./data/CDD/tax_data_fixed.rds')
  
rm(superceded, superced_ids, superceded_sum, superced_tax)

### end of superceded fixing ----


## Fix Updated  sequences ----
# 1034665437
# 1632804275

replaceable_ids <- 
  summaries |> 
  filter(!is.na(replacedby)) |> 
  select(nuc_id, accessionversion, replacedby) |> 
  mutate(
    token = map(replacedby, ~{
      entrez_search(db = 'nuccore', term = .x, 
                    use_history = T)
    }),
    fasta = map(
      token,
      ~entrez_fetch(db = 'nuccore', 
                    web_history = .x$web_history, 
                    rettype = 'fasta')
    ),
    nuc_summary = map(
      token,
      ~entrez_summary(db = 'nuccore',
                      web_history = .x$web_history)
    )
  ) |> 
  mutate(
    nuc_ss = map(fasta, fasta_to_DNAss),
    nuc_name = map(nuc_ss, names),
    nuc_seq = map(nuc_ss, paste)
  ) |> 
  select(nuc_id, replacedby, nuc_name, nuc_seq, nuc_summary)

print(replaceable_ids)

replacement_ids <- 
  replaceable_ids |> 
  unnest_wider(nuc_summary) |> 
  mutate(old_nuc_id = nuc_id,
         nuc_id = gi,
         tax_id = taxid) |> 
  select(-replacedby)
replacement_ids  

replaced_seq <- 
  replacement_ids |> 
  select(nuc_id, nuc_name, nuc_seq)

replaced_sum <- 
  replacement_ids |> 
  select(-nuc_name, -nuc_seq, -old_nuc_id) |> 
  mutate(across(everything(), as.character))

replacement_id_df <- 
  replacement_ids  |> 
  transmute(new_id = nuc_id,
            nuc_id = c('1632804275', '1034665437')) |> 
  left_join(id_data) |> 
  transmute(nuc_id = as.character(new_id),
            cdd_id, prot_id, tax_id)


# fix ids and summaries datasets
summaries <- summaries |> 
  filter(!nuc_id %in% replaceable_ids$nuc_id) |> 
  bind_rows(replaced_sum)

id_data <- id_data |> 
  filter(!nuc_id %in% replaceable_ids$nuc_id) |> 
  bind_rows(replacement_id_df)

# join 4 new seqs together
new_seqs <- bind_rows(
  superceded_seq |> select(nuc_id, nuc_name, nuc_seq),
  replaced_seq |> 
    mutate(nuc_id = as.character(nuc_id))
  ) |> 
  unnest(cols = c(nuc_name, nuc_seq))

print(new_seqs)


rm(replaceable_ids, replacement_id_df, replacement_ids, replaced_seq, replaced_sum,  superceded_seq)



## Removed* seqs ----

skimr::skim(summaries)

# comments for removed seqs
summaries |> 
  filter(!is.na(comment)) |> 
  count(comment)

# how many of these are missing taxonomy?
summaries |> 
  filter(!is.na(comment)) |> 
  left_join(id_data |> select(nuc_id, tax_id)) |> 
  distinct() |> 
  select(contains('id'), comment) |> 
  pull(tax_id) |> 
  is.na() |> sum()

# most of the ones missing taxonomy are comment = 'removed' 
id_data |> 
  select(nuc_id, tax_id) |> 
  pull(tax_id) |> 
  is.na() |> sum()

## might as well remove since can't find what to replace with.
remove_these_nuc_ids <- 
  summaries |> 
  filter(!is.na(comment)) |> 
  pull(nuc_id)
  
# filtered id dataset
id_data_filtered <- 
  id_data |> 
  filter(!nuc_id %in% remove_these_nuc_ids)
  
map(id_data_filtered, ~sum(is.na(.x)))

# still 9 missing tax ids have them in the summary,
# to retrieve
new_taxonomy <- 
  id_data_filtered |> 
  filter(is.na(tax_id)) |> 
  left_join(summaries) |> 
  mutate(tax_id = taxid) |> 
  fetch_taxonomy(id = tax_id)

# append to tax_data file
read_rds('./data/CDD/tax_data_fixed.rds') |> 
  bind_rows(new_taxonomy) |> 
  distinct() |> 
  write_rds('./data/CDD/tax_data_fixed.rds')

# get replacement id rows for those 9 missing taxonomy
id_data_fix <- 
  id_data_filtered |> 
  filter(is.na(tax_id)) |> 
  left_join(summaries) |> 
  mutate(tax_id = taxid) |> 
  select(contains('_id'))

# remove and re-bind replace 9 rows that were 9 missing tax ids.
id_data_filtered_fixed <- 
  id_data_filtered |> 
  filter(!is.na(tax_id)) |> 
  bind_rows(id_data_fix)

# now none are missing the taxonomy!!!
id_data_filtered_fixed |> 
  filter(is.na(tax_id))
  
rm(id_data_filtered, id_data_fix, new_taxonomy, remove_these_nuc_ids)

### end of removed seqs ----

## Remove one remaining RNAs

id_data_filtered_fixed_no_rna <- 
  id_data_filtered_fixed |> 
  left_join(summaries |> select(nuc_id, moltype)) |> 
  filter(!moltype == 'rna') |> 
  select(-moltype)


write_rds(id_data_filtered_fixed_no_rna,
          './data/CDD/ids_filtered.rds')

# filter summaries to remove 'removed' seqs
summaries <- summaries |> 
  filter(nuc_id %in% id_data_filtered_fixed_no_rna$nuc_id)

write_rds(summaries, './data/CDD/nuc_summary_filtered.rds')



summaries |> skimr::skim()
