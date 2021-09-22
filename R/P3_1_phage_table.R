### Part 3: Phages ------------------------------------------------------------

#' This script does the parsing of the genome accession tables and their metadata for phages of bacteria and archaea
#' Then it downloads all the Genbank files in 9 batches
#' Then the genbank XML text is parsed into dataframes and joined
#' Then the models are used on the protein sequences from the CDS info to identify integrases. 

## Parse lists of Genomes --------------------------------------------------

# I got the virus genome listed in viral genome tables from NCBI at:
# https://www.ncbi.nlm.nih.gov/genomes/GenomesGroup.cgi?taxid=10239&host=bacteria
# and 
# https://www.ncbi.nlm.nih.gov/genomes/GenomesGroup.cgi?taxid=10239&host=archaea

library(tidyverse)
library(janitor)
library(lubridate)
library(rentrez)
library(progress)

# parse tables for bacterial and archaeal virus genome records
all_entries <- 
  bind_rows(
    read_tsv('data/ncbi_viral_genomes/ncbi_bacteria_viruses.tbl', skip = 1),
    read_tsv('data/ncbi_viral_genomes/ncbi_archaea_viruses.tbl', skip = 1)
  ) |> 
  clean_names() |> 
  mutate(genome_length = str_extract(genome_length, '[0-9]+') |> as.numeric(),
         genome_neighbors = str_extract(genome_length, '[0-9]+') |> as.numeric(),
         date_completed = mdy(date_completed),
         date_updated = mdy(date_updated),
         number_of_segments =  as.numeric(number_of_segments)
  ) |> 
  mutate(across(where(is.character), as_factor)) |> 
  distinct()

glimpse(all_entries)


# 4260 genomes with only a single accession
single_segments <- 
  all_entries |> 
  filter(!is.na(accession)) |> 
  filter(accession != '-')

glimpse(single_segments)


# 57 genomes have multiple segments under different accessions
# needs to be parsed differently and joined with their header row (with no accession)
multiple_segments <- 
  all_entries |> 
  filter(is.na(accession)) |> 
  select(genome) |> 
  separate(genome, sep = '   ',
           into = c('genome', 'genome_length', 'accession',
                    'number_of_proteins', 'genome_neighbors')
           ) |> 
  mutate(accession = str_trim(accession)) |> 
  mutate(across(
    .cols = c(genome_length, number_of_proteins, genome_neighbors),
    .fns = ~str_extract(.x, '[0-9]+') |> as.numeric()
    )
  ) 
glimpse(multiple_segments)

# keep the metadata for viruses with multiple genome segments from the header rows which have a '-' symbol in their accession field
multiple_segments <- 
  all_entries |> 
  filter(accession == '-') |> 
  select(-c(accession, genome_length, number_of_proteins, genome_neighbors)) |> 
  full_join(multiple_segments, by = 'genome')
glimpse(multiple_segments)


# join the single and multiple segment genomes to get list of 4317 genomes
virus_tbl <- 
  full_join(single_segments, multiple_segments) |> 
  mutate(number_of_segments = ifelse(is.na(number_of_segments), 1, number_of_segments))

glimpse(virus_tbl)  

# notes:
# 109 'incomplete' genomes
# not sure what 'genome_neighbors' data is
summary(virus_tbl)

rm(all_entries, single_segments, multiple_segments)


## Download GenBank records from NCBI----

library(xml2)

## progress bar for n downloads
make_pb <- function(n){
  progress_bar$new(
    format = "Downloading Genbank records: [:bar] :percent eta: :eta",
    total = n,
    clear = FALSE,
  )
}

# entrez queries search & download xml formatted genbank file
fetch_viral_genome <- function(accession){
  search <- rentrez::entrez_search(db =  'nuccore', term = accession, use_history = T)
  genbank_xml <- rentrez::entrez_fetch(
    db = 'nuccore', 
    web_history = search$web_history, 
    retmode = 'xml',
    parsed = T,
    api_key = '889fdb9786a14019a0a1257196a09ba4ba08',
    rettype = 'gb')
  return(genbank_xml)
}


#### Execute downloads -------------------------------------------------------

# 1st group of genomes
acc <- virus_tbl |> 
  select(genome, accession) |> 
  slice(1:500) 
pb <- make_pb(nrow(acc))
xml1 <- acc |> 
  mutate(xml = map(accession, ~{
    print(glue::glue('\nAccession: {.x}'))
    xml <- fetch_viral_genome(.x)
    pb$tick()
    return(xml)
  }))
write_rds(xml1, 'data/ncbi_viral_genomes/genbank_xml/xml1.rds')
beepr::beep()

# 2nd 
acc <- virus_tbl |> 
  select(genome, accession) |> 
  slice(501:1000) 
pb <- make_pb(nrow(acc))
xml2 <- acc |> 
  mutate(xml = map(accession, ~{
    print(glue::glue('\nAccession: {.x}'))
    xml <- fetch_viral_genome(.x)
    pb$tick()
    return(xml)
  }))
write_rds(xml2, 'data/ncbi_viral_genomes/genbank_xml/xml2.rds')
beepr::beep()

# 3rd 
acc <- virus_tbl |> 
  select(genome, accession) |> 
  slice(1001:1500) 
pb <- make_pb(nrow(acc))
xml3 <- acc |> 
  mutate(xml = map(accession, ~{
    print(glue::glue('\nAccession: {.x}'))
    xml <- fetch_viral_genome(.x)
    pb$tick()
    return(xml)
  }))
write_rds(xml3, 'data/ncbi_viral_genomes/genbank_xml/xml3.rds')
beepr::beep()

# 4th
acc <- virus_tbl |> 
  select(genome, accession) |> 
  slice(1501:2000) 
pb <- make_pb(nrow(acc))
xml4 <- acc |> 
  mutate(xml = map(accession, ~{
    print(glue::glue('\nAccession: {.x}'))
    xml <- fetch_viral_genome(.x)
    pb$tick()
    return(xml)
  }))
write_rds(xml4, 'data/ncbi_viral_genomes/genbank_xml/xml4.rds')
beepr::beep()

# 5th
acc <- virus_tbl |> 
  select(genome, accession) |> 
  slice(2001:2500) 
pb <- make_pb(nrow(acc))
xml5 <- acc |> 
  mutate(xml = map(accession, ~{
    print(glue::glue('\nAccession: {.x}'))
    xml <- fetch_viral_genome(.x)
    pb$tick()
    return(xml)
  }))
write_rds(xml5, 'data/ncbi_viral_genomes/genbank_xml/xml5.rds')
beepr::beep()

# 6th
acc <- virus_tbl |> 
  select(genome, accession) |> 
  slice(2501:3000) 
pb <- make_pb(nrow(acc))
xml6 <- acc |> 
  mutate(xml = map(accession, ~{
    print(glue::glue('\nAccession: {.x}'))
    xml <- fetch_viral_genome(.x)
    pb$tick()
    return(xml)
  }))
write_rds(xml6, 'data/ncbi_viral_genomes/genbank_xml/xml6.rds')
beepr::beep()

# 7th
acc <- virus_tbl |> select(genome, accession) |> slice(3001:3500) 
pb <- make_pb(nrow(acc))
xml7 <- acc |> 
  mutate(xml = map(accession, ~{
    print(glue::glue('\nAccession: {.x}'))
    xml <- fetch_viral_genome(.x)
    pb$tick()
    return(xml)
  }))
write_rds(xml7, 'data/ncbi_viral_genomes/genbank_xml/xml7.rds')
beepr::beep()


# 8th
acc <- virus_tbl |> select(genome, accession) |> slice(3501:4000) 
pb <- make_pb(nrow(acc))
xml8 <- acc |> 
  mutate(xml = map(accession, ~{
    print(glue::glue('\nAccession: {.x}'))
    xml <- fetch_viral_genome(.x)
    pb$tick()
    return(xml)
  }))
write_rds(xml8, 'data/ncbi_viral_genomes/genbank_xml/xml8.rds')
beepr::beep()



# 9th
acc <- virus_tbl |> select(genome, accession) |> slice(4000:nrow(virus_tbl)) 
pb <- make_pb(nrow(acc))
xml9 <- acc |> 
  mutate(xml = map(accession, ~{
    print(glue::glue('\nAccession: {.x}'))
    xml <- fetch_viral_genome(.x)
    pb$tick()
    return(xml)
  }
  ))
write_rds(xml9, 'data/ncbi_viral_genomes/genbank_xml/xml9.rds')
beepr::beep()



## Parse Genbank Files -------------------------------------------------------

gb <- xmls$xml[[2]]

gbk <- 
  # xml to tibble
  read_xml(gb) |>
  as_list() |> 
  as_tibble() |> 
  # start flattening lists with only a single field (not features table)
  unnest_wider(GBSet) |> 
  rename_with(~str_remove(.x, 'GBSeq_'), .cols = everything()) |> 
  clean_names() |> 
  # ditch this garbage
  select(-c(xrefs, references, comment, keywords, other_seqids,
            matches('date'), matches('accession'))) |> 
  # flatten everything except for the features
  unnest(-feature_table) |> unnest(-feature_table)
  # mutate(feature_table = enframe(feature_table)) 
  
# get to features table and parse into nested table
ft <- gbk |> 
  select(locus, feature_table) |> 
  unnest_longer(feature_table, indices_include = F) |> 
  hoist(feature_table, 'GBFeature_key', .simplify = T) |> 
  clean_names() |> 
  unnest_longer(gb_feature_key) 

# go into feature table and parse the CDS features only
cds <- 
  ft |> 
  filter(gb_feature_key == 'CDS') |> 
  mutate(feat_id = row_number()) |> 
  hoist(feature_table, 'GBFeature_intervals') |>
  hoist(GBFeature_intervals, 'GBInterval') |> 
  hoist(GBInterval, 'GBInterval_from') |> 
  hoist(GBInterval, 'GBInterval_to') |> 
  unnest(c(GBInterval_from, GBInterval_to)) |> 
  select(-GBInterval) |> 
  unnest_longer(gb_feature_key) |> 
  hoist(feature_table, 'GBFeature_quals') |> 
  unnest(GBFeature_quals) |> 
  hoist(GBFeature_quals, 'GBQualifier_name') |>
  hoist(GBFeature_quals, 'GBQualifier_value') |> 
  unnest(c(GBQualifier_name, GBQualifier_value)) |> 
  select(-feature_table, contains('function')) |> 
  pivot_wider(names_from = GBQualifier_name, 
              values_from = GBQualifier_value,
              values_fn = list) |> 
  clean_names() |> 
  select(
    -contains('ribosomal_slippage'),
    -contains('inference'),
    -contains('function'),
    -contains('transl_except'),
    -contains('old_locus_tag'),
    -contains('note'), 
    -contains('ec_number'), 
    -contains('gene'), 
    -contains('transl_except'),
    -contains('artificial_location')
  )

ft_unnest <- cds |> 
  unnest(names(cds)[2:ncol(cds)]) |> 
  select(
    -contains('old_locus_tag'),
    -contains('note'), 
    -contains('ec_number'), 
    -contains('gene'), 
    -contains('transl_except'),
    -contains('artificial_location'),
    -contains('pseudo')
  )

ft_nest <- ft_unnest |> 
  unnest(gb_feature_key:translation) |> 
  nest(feature_table = gb_feature_key:translation)
tbl <- gbk |> 
  select(-feature_table) |> 
  distinct()
gb_record <- left_join(tbl, ft_nest, by = "locus")

# no CDS features in record, return an NA
if (nrow(ft) == 0) return(NA)


View(ft)
