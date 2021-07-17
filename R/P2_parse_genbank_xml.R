library(tidyverse)
library(rentrez)
library(xml2)
library(janitor)


# get genbank record as xml, parse xml
fetch_genbank <- function(nuc_id){
  # get record from entrez
  es <- entrez_search(db = 'nuccore', term = nuc_id, use_history = T)
  gbk <- entrez_fetch(web_history = es$web_history,
                      db = 'nuccore', rettype = 'gb', retmode = 'xml')
}

parse_genbank <- function(gbk){
  gbk <- read_xml(gbk) |> as_list()
  # start extraction
  tbl <- 
    tibble(GBSet=gbk) |> 
    unnest_wider(GBSet) |> 
    unnest_wider(GBSeq) |> 
    clean_names() |> 
    rename_all(~str_remove(.x, 'gb_seq_'))
  # start unnesting features
  tbl_unnest <- tbl |> 
    select(-c(contains('references'), 
              contains('keywords'),
              contains('other_seqids'),
              contains('contig'),
              contains('xrefs'),
              contains('sequence')
              )
           ) |> 
    relocate(-feature_table) |> 
    unnest(-feature_table) |> 
    unnest(-feature_table) |> 
    distinct()
  # feature table CDS extraction
  ft <- tbl_unnest |> 
    select(locus, feature_table) |> 
    unnest_longer(feature_table, indices_include = F) |> 
    hoist(feature_table, 'GBFeature_key', .simplify = T) |> 
    clean_names() |> 
    unnest_longer(gb_feature_key) |> 
    filter(gb_feature_key == 'CDS') 
  
  # no CDS features in record
  if (nrow(ft) == 0) return(NA)
  # full features extraction for CDS items in table only
  ft <- ft |> 
    mutate(feat_id = row_number()) |> 
    hoist(feature_table, 'GBFeature_intervals') |>
    hoist(GBFeature_intervals, 'GBInterval') |> 
    hoist(GBInterval, 'GBInterval_from') |> 
    hoist(GBInterval, 'GBInterval_to') |> 
    unnest(c(GBInterval_from, GBInterval_to)) |> 
    select(-GBInterval) |> 
    unnest_longer(gb_feature_key) |> 
    relocate(feature_table) |> 
    hoist(feature_table, 'GBFeature_quals') |> 
    unnest(GBFeature_quals) |> 
    hoist(GBFeature_quals, 'GBQualifier_name') |>
    hoist(GBFeature_quals, 'GBQualifier_value') |> 
    unnest(c(GBQualifier_name, GBQualifier_value)) |> 
    select(-feature_table) |> 
    pivot_wider(names_from = GBQualifier_name, values_from = GBQualifier_value) |> 
    clean_names() 
  ft <- ft |> 
    select(-contains('pseudo'),
           -contains('ribosomal_slippage'),
           -contains('inference'),
           -contains('function'),
           -contains('transl_except'),
           -contains('artificial_location'),
           -contains('old_locus_tag'),
           -contains('note'), 
           -contains('ec_number'), 
           -contains('gene'), 
           -contains('transl_except'),
           -contains('artificial_location'),
           -contains('pseudo')
           )
  ft_unnest <- ft |> 
    unnest(names(ft)[2:ncol(ft)]) |> 
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
  tbl <- tbl_unnest |> 
    select(-feature_table) |> 
    distinct()
  gb_record <- left_join(tbl, ft_nest, by = "locus")
  return(gb_record)
}




parse_genbank2 <- function(gbk){
  gbk <- read_xml(gbk) |> as_list()
  # start extraction
  tbl <- 
    tibble(GBSet=gbk) |> 
    unnest_wider(GBSet) |> 
    unnest_wider(`Bioseq-set_seq-set`) |>
    unnest_wider(`Seq-entry`) |>
    unnest_wider(`Seq-entry_seq`) |> 
    unnest_wider(`Bioseq`) |> View()
    unnest_wider(`Bioseq_id`)
    clean_names() |> 
    rename_all(~str_remove(.x, 'gb_seq_'))
  # start unnesting features
  tbl_unnest <- tbl |> 
    select(-c(contains('references'), 
              contains('keywords'),
              contains('other_seqids'),
              contains('contig'),
              contains('xrefs'),
              contains('sequence')
    )
    ) |> 
    relocate(-feature_table) |> 
    unnest(-feature_table) |> 
    unnest(-feature_table) |> 
    distinct()
  # feature table CDS extraction
  ft <- tbl_unnest |> 
    select(locus, feature_table) |> 
    unnest_longer(feature_table, indices_include = F) |> 
    hoist(feature_table, 'GBFeature_key', .simplify = T) |> 
    clean_names() |> 
    unnest_longer(gb_feature_key) |> 
    filter(gb_feature_key == 'CDS') 
  
  # no CDS features in record
  if (nrow(ft) == 0) return(NA)
  # full features extraction for CDS items in table only
  ft <- ft |> 
    mutate(feat_id = row_number()) |> 
    hoist(feature_table, 'GBFeature_intervals') |>
    hoist(GBFeature_intervals, 'GBInterval') |> 
    hoist(GBInterval, 'GBInterval_from') |> 
    hoist(GBInterval, 'GBInterval_to') |> 
    unnest(c(GBInterval_from, GBInterval_to)) |> 
    select(-GBInterval) |> 
    unnest_longer(gb_feature_key) |> 
    relocate(feature_table) |> 
    hoist(feature_table, 'GBFeature_quals') |> 
    unnest(GBFeature_quals) |> 
    hoist(GBFeature_quals, 'GBQualifier_name') |>
    hoist(GBFeature_quals, 'GBQualifier_value') |> 
    unnest(c(GBQualifier_name, GBQualifier_value)) |> 
    select(-feature_table) |> 
    pivot_wider(names_from = GBQualifier_name, values_from = GBQualifier_value) |> 
    clean_names() 
  ft <- ft |> 
    select(-contains('pseudo'),
           -contains('ribosomal_slippage'),
           -contains('inference'),
           -contains('function'),
           -contains('transl_except'),
           -contains('artificial_location'),
           -contains('old_locus_tag'),
           -contains('note'), 
           -contains('ec_number'), 
           -contains('gene'), 
           -contains('transl_except'),
           -contains('artificial_location'),
           -contains('pseudo')
    )
  ft_unnest <- ft |> 
    unnest(names(ft)[2:ncol(ft)]) |> 
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
  tbl <- tbl_unnest |> 
    select(-feature_table) |> 
    distinct()
  gb_record <- left_join(tbl, ft_nest, by = "locus")
  return(gb_record)
}



# id_data <- read_rds('./data/CDD/id_data_fixed.rds')
# nuc_summary <- read_rds('./data/CDD/nuc_summary_fixed.rds')
# # 
# # nuc_ids <- id_data |> 
# #   select(nuc_id) |> 
# #   distinct()
# # 
# # n <- nrow(nuc_ids)
# # chunk <- 100
# # map(seq(1, ))
# 
# # select one id to get genbank id
# nuc_id <- nuc_summary$nuc_id[19]
# gbk <- fetch_genbank(nuc_id)
# 
# nuc_summary[1002,1]
# # 
# gb <- fetch_genbank(nuc_id)
# # 
# gb
# gb |> unnest(feature_table)
# gb |> unnest(feature_table) |> View()

# 
# # get genbank recd
# es <- entrez_search(db = 'nuccore', term = s1, use_history = T)
# gbf <- entrez_fetch(web_history = es$web_history,
#                     db = 'nuccore', rettype = 'gb', retmode = 'xml')
# 
# ## don't need to write
# # gbk
# # tf <- tempfile()
# # write(gbk, tf)
# # x <- read_xml(gbk)
# # str(x)
# # x
# # x |>  xml_find_all('GBSet')        
# # str(x)
# # X <- as_list(x)
# # X
# 
# gbk <- read_xml(gbf) |> as_list()
# tbl <- 
#   tibble(GBSet=gbk) |> 
#   unnest_wider(GBSet) |> 
#   unnest_wider(GBSeq) |> 
#   clean_names() |> 
#   rename_all(~str_remove(.x, 'gb_seq_'))
# 
# tbl_unnest <- tbl |> 
#   select(-c(references, keywords, other_seqids, xrefs)) |> 
#   relocate(-feature_table) |> 
#   unnest(-feature_table) |> 
#   unnest(-feature_table) |> 
#   distinct()
# 
# ## THIS WORKS
# ft <- tbl_unnest |> 
#   select(locus, feature_table) |> 
#   unnest_longer(feature_table,indices_include = F) |> 
#   hoist(feature_table, 'GBFeature_key', .simplify = T) |> 
#   clean_names() |> 
#   unnest_longer(gb_feature_key) |> 
#   filter(gb_feature_key == 'CDS') |> 
#   mutate(feat_id = row_number()) |> 
#   hoist(feature_table, 'GBFeature_intervals') |>
#   hoist(GBFeature_intervals, 'GBInterval') |> 
#   hoist(GBInterval, 'GBInterval_from') |> 
#   hoist(GBInterval, 'GBInterval_to') |> 
#   unnest(c(GBInterval_from, GBInterval_to)) |> 
#   select(-GBInterval) |> 
#   unnest_longer(gb_feature_key) |> 
#   relocate(feature_table) |> 
#   hoist(feature_table, 'GBFeature_quals') |> 
#   unnest(GBFeature_quals) |> 
#   hoist(GBFeature_quals, 'GBQualifier_name') |>
#   hoist(GBFeature_quals, 'GBQualifier_value') |> 
#   unnest(c(GBQualifier_name, GBQualifier_value)) |> 
#   select(-feature_table) |> 
#   pivot_wider(names_from = GBQualifier_name, values_from = GBQualifier_value) |> 
#   clean_names() |> 
#   unnest(gb_feature_key:translation)
# 
# ft_nest <- ft |> 
#   unnest(gb_feature_key:translation) |> 
#   nest(feature_table = gb_feature_key:translation)
# 
# ft_nest |> unnest(col = feature_table) 
# 
# gb_record <- left_join(tbl_unnest |> select(-feature_table), ft_nest)
# gb_record 
# gb_record |> unnest(col = feature_table) |> View()
# 
# 
# 
# 
# 
# 
# 
