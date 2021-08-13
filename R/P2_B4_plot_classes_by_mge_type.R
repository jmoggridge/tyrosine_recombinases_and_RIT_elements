## Nicole wants to know if certain integrase types are more common in IME
library(tidyverse)
library(glue)
# devtools::install_github("hrbrmstr/waffle")
library(waffle)
library(rcartocolor)

fix_this <- 
  read_rds('./data/iceberg/iceberg_db_classed_proteins.rds') |> 
  ungroup() |> 
  filter(!is.na(prot_seq)) |> 
  # remove these improperly identified ICEs
  filter(str_detect(name, 'ICEMlSym')) |> 
  mutate(across(c(knn, rf, glmnet), ~as.character(.x))) |> 
  mutate(
    mge_name = name,
    consensus_pred = case_when(
      knn == glmnet ~ knn,
      knn == rf ~  knn,
      rf == glmnet ~ rf,
      rf != knn & rf != glmnet ~'no consensus',
      TRUE ~ 'NA'
    ),
    description = str_remove(description, '^.*?\\s')
  ) |> 
  select(-name, -dna_seq) |> 
  relocate(type, id, mge_name, organism) |> 
  distinct() |> 
  filter(consensus_pred != 'Other')

iceberg <- 
  read_rds('./data/iceberg/iceberg_db_classed_proteins.rds') |> 
  ungroup() |> 
  filter(!is.na(prot_seq)) |> 
  # remove these improperly identified ICEs
  filter(!name == 'ICEMlSym(NZP2037)') |> 
  # get consensus of predicted classes
  mutate(across(c(knn, rf, glmnet), ~as.character(.x))) |> 
  mutate(
    mge_name = name,
    consensus_pred = case_when(
      knn == glmnet ~ knn,
      knn == rf ~  knn,
      rf == glmnet ~ rf,
      rf != knn & rf != glmnet ~'no consensus',
      TRUE ~ 'NA'
    ),
    description = str_remove(description, '^.*?\\s')
  ) |> 
  select(-name) |> 
  relocate(type, id, mge_name, organism)

glimpse(iceberg)

# number of elements by types
mge_count <- iceberg |> 
  select(type, mge_name, id) |> 
  distinct() |> 
  dplyr::count(type) 
mge_count
  
# number of proteins by mge types
mge_proteins_count <- iceberg |> 
  count(type, name = 'n_proteins')
mge_proteins_count

# total number of integrases across types
mge_integrase_count <- iceberg |> 
  filter(!consensus_pred %in% c('no consensus', 'Other')) |> 
  count(type, name = 'n_integrases')
mge_integrase_count



# first bar plot for Nicole: ICEberg integrases counts by mge type, without mislabelled tripartite ice
iceberg_mge_integrases_barplot <- iceberg |> 
  left_join(mge_integrase_count) |> 
  mutate(type = glue::glue('{type} ({n_integrases} integrases)')) |> 
  filter(!consensus_pred %in% c('no consensus', 'Other')) |> 
  group_by(type) |> 
  ggplot(aes(y = fct_rev(consensus_pred))) + 
  geom_bar(show.legend = F) +
  facet_wrap(~type) +
  theme_bw() +
  labs(y = 'Predicted integrase subfamily', 
       fill = 'MGE type\n(n integrases)',
       title = 'Integrase subfamilies across MGE types',
       subtitle = 'ICEberg2.0 - experimental, intact database')
iceberg_mge_integrases_barplot

write_rds(iceberg_mge_integrases_barplot, './results/iceberg_mge_integrases_barplot.rds')

## alternately with stacked bars



## for questions pertaining to number of integrases...

# nest all proteins data for each MGE
iceberg_nested <- iceberg |> 
  select(-c(parent_element, accession, organism, description, start_pos, end_pos,
            matches('seq|mystery|parent|database|prob'))) |> 
  group_by(type, id, mge_name) |> 
  nest(prot_data = c(matches('prot_'), knn, rf, glmnet, consensus_pred)) |> 
  ungroup()
glimpse(iceberg_nested)


iceberg_integrase_counts <- 
  iceberg_nested |> 
  mutate(integrase_count = map(
    prot_data, 
    ~{.x |> 
        count(consensus_pred) |> 
        mutate(consensus_pred = glue('count_{consensus_pred}')) |> 
        pivot_wider(names_from = consensus_pred, values_from = n)
    })) |>
  unnest_wider(integrase_count) |> 
  select(-count_Other, -`count_no consensus`) |> 
  mutate(across(starts_with('count_'), ~ifelse(is.na(.x), 0, .x))) |>
  mutate(n_integrases = rowSums(across(starts_with('count_'))),
         group = case_when(
           n_integrases == 0 ~ 'no integrase',
           n_integrases == 1 ~ '1 integrase',
           n_integrases > 1 ~ '2+ integrases',
           TRUE ~ 'NA'
         ))

glimpse(iceberg_integrase_counts)
# iceberg_integrase_counts |> View()
# iceberg_integrase_counts |> filter(type == 'CIME') |> View()

## How many integrases are present per element?
iceberg_integrase_counts |> 
  filter(!type == 'T4SS-type_ICE') |> 
  ggplot(aes(x = n_integrases)) +
  geom_bar() + 
  theme_bw() +
  facet_wrap(~type, nrow = 1, scales = 'free_y') +
  labs(x = 'Number of integrases per MGE', y = 'Number of MGEs')

# what integrases are present on elements with only a single integrase/ 2+ integrases?
counts_by_element <- 
  iceberg_integrase_counts |> 
  filter(group != 'no integrase') |> 
  pivot_longer(cols = count_Int_Tn916:count_RitC, 
               names_to = 'subfamily', values_to = 'count') |> 
  mutate(subfamily = fct_rev(str_remove(subfamily, 'count_')))

counts_by_element |> 
  filter(!type== 'T4SS-type_ICE') |> 
  ggplot(aes(y = subfamily, x = count)) +
  geom_col() +
  facet_grid(group~type) +
  theme_bw() +
  labs(x = 'Count', y = 'Integrase subfamily')


# plot of co-occurrence of subfamilies on MGEs for ICEs and IMEs.
counts_by_element |> 
  filter(group == '2+ integrases') |> 
  filter(!type == 'T4SS-type_ICE') |> 
  ggplot(aes(y = fct_rev(mge_name), x = count, fill = fct_rev(subfamily))) +
  geom_col(position = 'stack') +
  scale_fill_carto_d() +
  facet_grid(type~., scales = 'free_y', space = 'free_y') +
  labs(y = 'MGE name', x = 'Count', fill = 'subfamily') +
  theme_bw()
  

# Do we get ICEs with both a Tn916 and an IntSXT? NO!
# I would expect the integrons and RIT elements to only occur on elements that also carry one of the putative functional groups since they only mobilize internally within a cell and wouldn't be responsible for moving the ICE. 


