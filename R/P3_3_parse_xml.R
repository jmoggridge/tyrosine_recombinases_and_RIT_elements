## Part 3, step 2: Viral genomes 
library(tidyverse)
library(furrr)
library(janitor)
library(xml2)
library(progress)

# genbank parsing function adapted from earlier script P2_parse_genbank.R
parse_genbank_xml <- function(gb){
  gbk <- 
    # xml to tibble
    read_xml(gb) |>
    as_list() |> 
    as_tibble() |> 
    # start flattening lists with only a single field (not features table)
    unnest_wider(GBSet) |> 
    rename_with(~str_remove(.x, 'GBSeq_'), .cols = everything()) |> 
    janitor::clean_names() |> 
    # ditch this garbage
    select(-c(contains('xrefs'), 
              contains('references'),
              contains('comment'),
              contains('keywords'),
              contains('other_seqids'),
              contains('experiment'),
              contains('citation'),
              matches('date'), 
              matches('accession')
    )) |> 
    # flatten everything except for the features
    unnest(-feature_table) |> unnest(-feature_table)
  # mutate(feature_table = enframe(feature_table)) 
  
  # get to features table and parse into nested table
  ft <- gbk |> 
    select(locus, feature_table) |> 
    unnest_longer(feature_table, indices_include = F) |> 
    hoist(feature_table, 'GBFeature_key', .simplify = T) |> 
    janitor::clean_names() |> 
    unnest_longer(gb_feature_key) 

  # go into feature table and parse the CDS features only
  cds <- ft |> 
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
    janitor::clean_names() |> 
    select(
      -contains('ribosomal_slippage'),
      -contains('inference'),
      -contains('function'),
      -contains('transl_except'),
      -contains('old_locus_tag'),
      -contains('db_xref'),
      -contains('locus_tag'),
      -contains('note'), 
      -contains('ec_number'), 
      -contains('gene'), 
      -contains('transl_except'),
      -contains('artificial_location')
    )
  
  # no CDS features in record, return an NA
  if (nrow(ft) == 0) return(NA)
  
  ft_unnest <- cds |> 
    unnest(names(cds)[2:ncol(cds)]) |> 
    select(
      -contains('old_locus_tag'),
      -contains('note'), 
      -contains('ec_number'), 
      -contains('gene'), 
      -contains('transl_except'),
      -contains('artificial_location'),
      -contains('experiment'),
      -contains('citation'),
      -contains('exception'),
      -contains('standard_name'),
      -contains('number'),
      -contains('pseudo')
    )
  
  ft_nest <- ft_unnest |> 
    unnest(gb_feature_key:translation) |> 
    nest(feature_table = gb_feature_key:translation)
  tbl <- gbk |> 
    select(-feature_table) |> 
    distinct()
  gb_record <- left_join(tbl, ft_nest, by = "locus")
  
  return(gb_record)
}

# parse genbank w printed accession for each step
parse_gb_xml_verbose <- function(gb, accession){
  cat(glue::glue('{accession}\n'))
  parse_genbank_xml(gb)
}


## Load raw genbank records -------------------------------------------------

# load previously downloaded genbank xml records.
virus_df <- 
  Sys.glob('data/ncbi_viral_genomes/genbank_xml/*') |> 
  map(read_rds) |> 
  reduce(bind_rows) |> 
  # get rid of duplicates
  distinct()

glimpse(virus_df)

# split up un-parsed genbank entries 
virus1 <- virus_df |> slice(1:1000)
virus2 <- virus_df |> slice(1001:2000)
virus3 <- virus_df |> slice(2001:3000)
virus4 <- virus_df |> slice(3001:4000)
virus5 <- virus_df |> slice(4001:nrow(virus_df))


## Parse Genbank recds in batches
plan(multisession, workers = 8)
gb <- virus_df |>
  mutate(parsed = future_map2(
    .x = xml, .y = accession, .f = parse_gb_xml_verbose, .progress = T
  )) |> 
  select(-xml) |> 
  unnest(parsed) |> 
  distinct() |> 
  select(-sequence, -contains('primary'))
beepr::beep()

gb <- gb |> select(-primary)
# gb1 |> select(-feature_table) |> View()
write_rds(gb, 'data/ncbi_viral_genomes/genbank_parsed.rds')

# 
# #1000
# plan(multisession, workers = 8)
# gb1 <- virus1 |>
#   mutate(parsed = future_map2(
#     .x = xml, .y = accession, .f = parse_gb_xml_verbose, .progress = T
#   )) |> 
#   select(-xml) |> 
#   unnest(parsed) |> 
#   distinct() |> 
#   select(-sequence)
# beepr::beep()
# 
# # gb1 |> select(-feature_table) |> View()
# write_rds(gb1, 'data/ncbi_viral_genomes/genbank_parsed/gb1.rds')
# rm(virus1)
# 
# #2000
# plan(multisession, workers = 8)
# gb2 <- virus2 |>
#   mutate(parsed = future_map2(
#     .x = xml, .y = accession, .f = parse_gb_xml_verbose, .progress = T
#   )) |> 
#   select(-xml) |> 
#   unnest(parsed) |> 
#   distinct() |> 
#   select(-sequence)
# beepr::beep()
# rm(virus2)
# write_rds(gb2, 'data/ncbi_viral_genomes/genbank_parsed/gb2.rds')
# 
# #3000
# gb3 <- virus3 |>
#   mutate(parsed = future_map2(
#     .x = xml, .y = accession, .f = parse_gb_xml_verbose, .progress = T
#   )) |> 
#   select(-xml) |> 
#   unnest(parsed) |> 
#   distinct() |> 
#   select(-sequence)
# beepr::beep()
# rm(virus3)
# write_rds(gb3, 'data/ncbi_viral_genomes/genbank_parsed/gb3.rds')
# 
# #4000
# gb4 <- virus4 |>
#   mutate(parsed = future_map2(
#     .x = xml, .y = accession, .f = parse_gb_xml_verbose, .progress = T
#   )) |> 
#   select(-xml) |> 
#   unnest(parsed) |> 
#   distinct() |> 
#   select(-sequence)
# beepr::beep()
# rm(virus4)
# write_rds(gb4, 'data/ncbi_viral_genomes/genbank_parsed/gb4.rds')
# 
# 
# #5000
# gb5 <- virus5 |>
#   mutate(parsed = future_map2(
#     .x = xml, .y = accession, .f = parse_gb_xml_verbose, .progress = T
#   )) |> 
#   select(-xml) |> 
#   unnest(parsed) |> 
#   distinct() |> 
#   select(-sequence)
# beepr::beep()
# rm(virus5)
# # gb5 |> select(-feature_table) |> View()
# 
# write_rds(gb5, 'data/ncbi_viral_genomes/genbank_parsed/gb5.rds')
# 
# 
# viruses <- bind_rows(gb1, gb2, gb3, gb4, gb5)
# 
# rm(gb1, gb2, gb3, gb4, gb5)



# ## ISSUES WITH THESE - WEIRD FIELDS I HAVEN'T SEEN BEFORE
# # has db_xref and citation columns
# gb <- virus1 |> 
#   filter(accession == 'NC_049976') |> 
#   pull(xml) |> 
#   pluck(1)
# parse_genbank_xml(gb)
# 
# gb <- virus1 |>
#   filter(accession == 'NC_004166') |>
#   pull(xml) |>
#   pluck(1)
# parse_genbank_xml(gb)

# gb <- virus2 |>
#   filter(accession == 'NC_054637') |>
#   pull(xml) |>
#   pluck(1)
# parse_genbank_xml(gb)

# rm(cds, ft, ft_nest, ft_unnest, gbk)
