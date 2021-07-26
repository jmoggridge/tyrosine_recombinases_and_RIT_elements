library(tidyverse)
library(rentrez)
library(xml2)
library(janitor)
library(progress)

## progress bar for n iters
make_pb <- function(n){
  progress_bar$new(
    format = "Finding Rits: [:bar] :percent eta: :eta",
    total = n,
    clear = FALSE,
  )
}

# get genbank record as xml
fetch_genbank <- function(nuc_id){
  # get record from entrez
  es <- entrez_search(db = 'nuccore', term = nuc_id, use_history = T)
  gbk <- entrez_fetch(web_history = es$web_history,
                      db = 'nuccore', rettype = 'gb', retmode = 'xml')
}

# get genbank WITH PARTS record as xml, parse xml
fetch_genbank_with_parts <- function(nuc_id){
  # get record from entrez
  es <- entrez_search(db = 'nuccore', term = nuc_id, use_history = T)
  gbk <- entrez_fetch(web_history = es$web_history,
                      db = 'nuccore', rettype = 'gbwithparts', retmode = 'xml')
}


# parse vanilla 'gb' type from mode 'xml'
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
    pivot_wider(names_from = GBQualifier_name, 
                values_from = GBQualifier_value,
                values_fn = list) |> 
    clean_names() 
  ft <- ft |> 
    select(
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
           -contains('artificial_location')
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

# 
# 
# parse_genbank_huge <- function(gbk){
#   gbk <- read_xml(gbk, options = 'HUGE', ) |> as_list()
#   # start extraction
#   tbl <- 
#     tibble(GBSet=gbk) |> 
#     unnest_wider(GBSet) |> 
#     select(GBSeq) |> 
#     unnest_wider(`GBSeq`) |>
#     unnest_wider(`Seq-entry`) |>
#     unnest_wider(`Seq-entry_seq`) |> 
#     unnest_wider(`Bioseq`)
#   unnest_wider(`Bioseq_id`)
#   clean_names() |> 
#     rename_all(~str_remove(.x, 'gb_seq_'))
#   # start unnesting features
#   tbl_unnest <- tbl |> 
#     select(-c(contains('references'), 
#               contains('keywords'),
#               contains('other_seqids'),
#               contains('contig'),
#               contains('xrefs'),
#               contains('sequence')
#     )
#     ) |> 
#     relocate(-feature_table) |> 
#     unnest(-feature_table) |> 
#     unnest(-feature_table) |> 
#     distinct()
#   # feature table CDS extraction
#   ft <- tbl_unnest |> 
#     select(locus, feature_table) |> 
#     unnest_longer(feature_table, indices_include = F) |> 
#     hoist(feature_table, 'GBFeature_key', .simplify = T) |> 
#     clean_names() |> 
#     unnest_longer(gb_feature_key) |> 
#     filter(gb_feature_key == 'CDS') 
#   
#   # no CDS features in record
#   if (nrow(ft) == 0) return(NA)
#   # full features extraction for CDS items in table only
#   ft <- ft |> 
#     mutate(feat_id = row_number()) |> 
#     hoist(feature_table, 'GBFeature_intervals') |>
#     hoist(GBFeature_intervals, 'GBInterval') |> 
#     hoist(GBInterval, 'GBInterval_from') |> 
#     hoist(GBInterval, 'GBInterval_to') |> 
#     unnest(c(GBInterval_from, GBInterval_to)) |> 
#     select(-GBInterval) |> 
#     unnest_longer(gb_feature_key) |> 
#     relocate(feature_table) |> 
#     hoist(feature_table, 'GBFeature_quals') |> 
#     unnest(GBFeature_quals) |> 
#     hoist(GBFeature_quals, 'GBQualifier_name') |>
#     hoist(GBFeature_quals, 'GBQualifier_value') |> 
#     unnest(c(GBQualifier_name, GBQualifier_value)) |> 
#     select(-feature_table) |> 
#     pivot_wider(names_from = GBQualifier_name, values_from = GBQualifier_value) |> 
#     clean_names() 
#   ft <- ft |> 
#     select(-contains('pseudo'),
#            -contains('ribosomal_slippage'),
#            -contains('inference'),
#            -contains('function'),
#            -contains('transl_except'),
#            -contains('artificial_location'),
#            -contains('old_locus_tag'),
#            -contains('note'), 
#            -contains('ec_number'), 
#            -contains('gene'), 
#            -contains('transl_except'),
#            -contains('artificial_location'),
#            -contains('pseudo')
#     )
#   ft_unnest <- ft |> 
#     unnest(names(ft)[2:ncol(ft)]) |> 
#     select(
#       -contains('old_locus_tag'),
#       -contains('note'), 
#       -contains('ec_number'), 
#       -contains('gene'), 
#       -contains('transl_except'),
#       -contains('artificial_location'),
#       -contains('pseudo')
#     )
#   
#   ft_nest <- ft_unnest |> 
#     unnest(gb_feature_key:translation) |> 
#     nest(feature_table = gb_feature_key:translation)
#   tbl <- tbl_unnest |> 
#     select(-feature_table) |> 
#     distinct()
#   gb_record <- left_join(tbl, ft_nest, by = "locus")
#   return(gb_record)
# }
