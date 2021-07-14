library(tidyverse)
library(rentrez)
library(jsonlite)
library(janitor)
library(xml2)

id_data <- read_rds('./data/CDD/id_data_fixed.rds')
nuc_summary <- read_rds('./data/CDD/nuc_summary_fixed.rds')

# select one id to get genbank id
s1 <- nuc_summary$nuc_id[1]
nuc_summary[1,1]
# get genbank recd
es <- entrez_search(db = 'nuccore', term = s1, use_history = T)
gbk <- entrez_fetch(web_history = es$web_history,
                    db = 'nuccore', rettype = 'gb', retmode = 'xml')

tf <- tempfile()
write(gbk, tf)

x <- read_xml(tf)
str(x)
x %>% xml_find_all('GBSet')        
str(x)
X <- as_list(x)
X
tbl <- 
  tibble(GBSet=X) |> 
  unnest_wider(GBSet) |> 
  unnest_wider(GBSeq) |> 
  clean_names() |> 
  rename_all(~str_remove(.x, 'gb_seq_'))

tbl_unnest <- tbl |> 
  select(-c(references, keywords, other_seqids, xrefs)) |> 
  relocate(-feature_table) |> 
  unnest(-feature_table) |> 
  unnest(-feature_table) |> 
  distinct()

## THIS WORKS
ft <- tbl_unnest |> 
  select(locus, feature_table) |> 
  unnest_longer(feature_table,indices_include = F) |> 
  hoist(feature_table, 'GBFeature_key', .simplify = T) |> 
  clean_names() |> 
  unnest_longer(gb_feature_key) |> 
  filter(gb_feature_key == 'CDS') |> 
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
  pivot_wider(names_from = GBQualifier_name, values_from = GBQualifier_value)

ft |> View()







