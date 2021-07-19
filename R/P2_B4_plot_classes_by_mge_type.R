## Nicole wants to know if certain integrase types are more common in IME
library(tidyverse)

iceberg <- 
  read_rds('./data/iceberg/iceberg_db_classed_proteins.rds') |> 
  filter(!is.na(prot_seq)) |> 
  select(-contains('seq'))
glimpse(iceberg)


# get consensus of predicted classes
iceberg <- iceberg |> 
  mutate(across(c(knn, rf, glmnet), ~as.character(.x))) |> 
  mutate(consensus_pred = case_when(
    knn == glmnet ~ knn,
    knn == rf ~  knn,
    rf == glmnet ~ rf,
    rf != knn & rf != glmnet ~'no consensus',
    TRUE ~ 'NA'
  )) 

# total number of elements across types
mge_count <- iceberg |> 
  select(type, name, id) |> 
  distinct() |> 
  count(type)
mge_count

# total number of proteins across types
mge_proteins_count <- iceberg |> 
  count(type, name = 'n_proteins')
mge_proteins_count

# total number of integrases across types
mge_integrase_count <- iceberg |> 
  filter(!consensus_pred %in% c('no consensus', 'Other')) |> 
  count(type, name = 'n_integrases')
mge_integrase_count

# bar plot for Nicole
iceberg_mge_integrases_barplot <- 
  iceberg |> 
  left_join(mge_integrase_count) |> 
  mutate(type = glue::glue('{type} ({n_integrases})')) |> 
  filter(!consensus_pred %in% c('no consensus', 'Other')) |> 
  group_by(type) |> 
  ggplot(aes(y = fct_rev(consensus_pred), fill = type)) + 
  geom_bar() +
  labs(y = 'Predicted integrase subfamily', 
       fill = 'MGE type\n(n integrases)',
       title = 'Integrase subfamilies across MGE types',
       subtitle = 'ICEberg2.0 - experimental, intact database')
  
write_rds(iceberg_mge_integrases_barplot, './results/iceberg_mge_integrases_barplot.rds')
