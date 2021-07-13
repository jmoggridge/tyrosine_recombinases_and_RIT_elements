## P2_A2_repair_nucleotide_data -----

library(tidyverse)
source('./R/P2_entrez_functions.R')

## read data ----

# ids dataset from P2_A1
id_data <- 
  read_rds('./data/CDD/cdd_prot_nuc_tax_ids.rds') |> 
  filter(!is.na(nuc_id))
glimpse(id_data)

# track missing ids and number of unique ids 
map(id_data |> select(nuc_id, tax_id) |> distinct(),
    ~sum(is.na(.x)))


summaries <- 
  read_rds('./data/CDD/nuc_summary.rds') |> 
  mutate(nuc_id = name) |> 
  select(-c(name, gi, uid, term)) |> # all same as nuc_id
  mutate(across(everything(), ~na_if(.x, '')))
# glimpse(summaries)
# skimr::skim(summaries)

## Issues? Some EDA ----

# join ids with their summaries
df <- id_data |>
  select(nuc_id) |> 
  distinct() |> 
  left_join(summaries)

# 207 comments
# 203 sequences removed, 2 superceded, 2 updated
df |> count(comment)
# 205 'suppressed'; 2 replaced
df |> count(status)
# updated: 1034665437, 1632804275
df |> 
  filter(!is.na(replacedby)) |> 
  count(nuc_id, replacedby, comment)

# most circular are probably complete genomes or maybe plasmids
df |> count(topology)
# 1 circular  1231
# 2 linear   10445
# 3 not-set      7

df |> count(genome)
# 1 chromosome   901
# 2 genomic     5041
# 3 plasmid      365
# 4 NA          5376

# two seqs are RNA... should exclude. possibly duplicate nuc records for same protein with RNA and DNA?
df |> count(moltype)
df |> count(biomol)

#  mostly wgs
df |> count(tech)
# ~ half 'refseq', half 'insd'
df |> count(sourcedb)

# # seq lengths seem okay
# df |>
#   transmute(nuc_id, kbp = as.numeric(slen)/1000) |> 
#   ggplot(aes(kbp)) +
#   geom_histogram() +
#   scale_x_log10() +
#   labs(title = glue('Lengths of distinct nucleotide seqs'))


## Superceded and updated replacement records -----
superceded_updated <- 
  bind_rows(
    # updated
    df |> 
      filter(!is.na(replacedby)) |> 
      mutate(new_nuc_id = replacedby),
    # superceded
    df |> 
      filter(str_detect(comment, 'CP053749|CP063356')) |> 
      mutate(new_nuc_id = str_extract(
        comment, 'CP053749|CP063356')
        )
  ) |> 
  left_join(id_data) |> 
  select(contains('_id')) |> 
  mutate(
    token = map(
      new_nuc_id, 
      ~entrez_search(db = 'nuccore', term = .x,
                     use_history = T)
    ),
    fasta = map(
      token,
      ~entrez_fetch(db = 'nuccore', 
                    web_history = .x$web_history, 
                    rettype = 'fasta')
    ),
    nuc_summary = map(
      token,
      ~ entrez_summary(db = 'nuccore',
                       web_history = .x$web_history)
    )
  ) |> 
  mutate(
    nuc_ss = map(fasta, fasta_to_DNAss),
    nuc_name = map(nuc_ss, names),
    nuc_seq = map(nuc_ss, paste)
  ) |> 
  select(contains('id'), contains('nuc_')) |> 
  unnest_wider(nuc_summary) |> 
  mutate(
    old_nuc_id = nuc_id,
    nuc_id = uid,
    tax_id = taxid
  ) |> 
  select(-c(gi, uid, term, new_nuc_id, nuc_ss)) |> 
  mutate(across(everything(), as.character))

superceded_updated |> 
  select(contains('_id'))

# remove old ids, replace with new rows
id_data_fix1 <- 
  id_data |> 
  filter(!nuc_id %in% superceded_updated$old_nuc_id) |> 
  bind_rows(superceded_updated |> 
              select(contains('_id')) |> 
              select(-old_nuc_id)
  )
glimpse(id_data_fix1)
# missing tax id down to 210 from 214 bc four fixed
map(id_data_fix1 |> select(nuc_id, tax_id) |> distinct(),
    ~sum(is.na(.x)))

# remove old summaries, replace with new rows
summaries_fix1 <- 
  summaries |> 
  filter(!nuc_id %in% superceded_updated$old_nuc_id) |> 
  bind_rows(superceded_updated) |> 
  select(-c(nuc_name, nuc_seq, old_nuc_id, tax_id,
            prot_id, cdd_id))

# replacement sequences
superceded_updated_seq <- superceded_updated |> 
  select(nuc_id, old_nuc_id, nuc_name, nuc_seq)

# replacement taxonomy
superceded_updated_tax <- superceded_updated |> 
  select(nuc_id, tax_id) |> 
  fetch_taxonomy(id = tax_id) |> 
  select(-nuc_id)
superceded_updated_tax

read_rds('./data/CDD/tax_data.rds') |> 
  bind_rows(superceded_updated_tax) |> 
  distinct() |> 
  write_rds('./data/CDD/tax_data_fixed.rds')
  
rm(superceded_updated, superceded_updated_tax, df,
   summaries, id_data)


## Suppressed records removal  ----
## just get rid of these entirely!

# comments for removed seqs
summaries_fix1 |> 
  filter(!is.na(comment)) |> 
  count(comment)
summaries_fix1 |> 
  filter(!is.na(comment)) |> 
  count(comment) |> 
  pull(comment)

# remove all 205 of those
to_remove <- 
  summaries_fix1 |> 
  filter(!is.na(comment))

id_data_fix2 <- 
  id_data_fix1 |> 
  filter(!nuc_id %in% to_remove$nuc_id)

summaries_fix2 <- 
  summaries_fix1 |> 
  filter(!nuc_id %in% to_remove$nuc_id)

rm(to_remove, id_data_fix1, summaries_fix1)

# now only a few ids are still missing taxonomy
map(id_data_fix2 |> 
      select(nuc_id, tax_id) |> 
      distinct(),
    ~sum(is.na(.x)))


# still 7 nuc id missing tax ids, but have their `taxid` in the summary, so to retrieve:
id_with_fixed_tax_id <- 
  id_data_fix2 |> 
  filter(is.na(tax_id)) |> 
  left_join(summaries_fix2, by = "nuc_id") |> 
  mutate(tax_id = taxid) |> 
  select(contains('_id'))

id_with_fixed_tax_id

# fetch missing tax records and append to tax_data file
fixed_taxonomy <- id_with_fixed_tax_id |> 
  fetch_taxonomy(id = tax_id) |> 
  select(-c(nuc_id, cdd_id, prot_id))
  
fixed_taxonomy

read_rds('./data/CDD/tax_data_fixed.rds') |> 
  bind_rows(fixed_taxonomy) |> 
  distinct() |> 
  write_rds('./data/CDD/tax_data_fixed.rds')

# create replacement id rows for those 9 missing taxonomy
id_data_fix3 <- 
  id_data_fix2 |> 
  filter(!nuc_id %in% id_with_fixed_tax_id$nuc_id) |> 
  bind_rows(id_with_fixed_tax_id)

# id_data is now fixed
id_data_fix3
map(id_data_fix3, ~sum(is.na(.x)))


# summaries for all nucs?
length(unique(summaries_fix2$nuc_id)) ==
  length(unique(id_data_fix3$nuc_id)) 

# yes ok, keep names matching
summaries_fix3 <- summaries_fix2

rm(id_with_fixed_tax_id, fixed_taxonomy,
   id_data_fix2, summaries_fix2)


## RNA record removal -----
is_rna <- 
  id_data_fix3 |> 
  left_join(summaries_fix3 |> select(nuc_id, moltype)) |> 
  filter(moltype == 'rna')

id_data_fix4 <- 
  id_data_fix3 |> 
  filter(!nuc_id %in% is_rna$nuc_id)

summaries_fix4 <- 
  summaries_fix3 |> 
  filter(!moltype == 'rna')

# now have summaries for all nucs?
all(
  sort(unique(summaries_fix4$nuc_id)) ==
      sort(unique(id_data_fix4$nuc_id)) 
  )
# still have complete set of ids for each nuc?
map(id_data_fix4, ~sum(is.na(.x)))

rm(id_data_fix3, summaries_fix3, is_rna)

# write fixed id_data and nuc_summaries
write_rds(id_data_fix4, './data/CDD/id_data_fixed.rds')
write_rds(summaries_fix4, './data/CDD/nuc_summary_fixed.rds')

# remove discarded records from prot_data
read_rds('./data/CDD/prot_data.rds') |> 
  filter(prot_id %in% id_data_fix4$prot_id) |> 
  write_rds('./data/CDD/prot_data_fixed.rds')

# removed discarded records from prot_summary 
read_rds('./data/CDD/prot_summary.rds') |> 
  select(-token) |> 
  mutate(summary = map(summary, ~{
    .x |> mutate(across(everything(), as.character))
  })) |> 
  unnest(c(id, summary)) |> 
  filter(id %in% id_data_fix4$prot_id) |> 
  mutate(prot_id = id) |> 
  select(-id) |> 
  write_rds('./data/CDD/prot_summary_fixed.rds')




