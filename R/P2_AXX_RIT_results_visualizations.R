## RIT distinct results EDA and grouping by similarity
library(tidyverse)

# full ungrouped data
rit_elements <- read_rds('./results/rit_elements.rds')
# grouped by tax_id & rit_dna sequence
distinct_rits <- read_rds('results/distinct_RITs.rds')


## EDA taxonomy ============================================================

# 187 distinct taxon ids; one taxon has 35! distinct RITs.
distinct_rits |> dplyr::count(tax_id) |> arrange(desc(n))

distinct_rits |>
  dplyr::count(tax_id) |> 
  ggplot(aes(n)) + 
  geom_histogram() + 
  labs(
    x = 'number of distinct elements',
    y = 'number of taxa',
    title = 'Distribution of the number of distinct elements per strain')

# one strain has 32 RITs. -> probably just SNPs??
distinct_rits |>
  dplyr::count(tax_name) |> 
  filter(n > 2) |> 
  ggplot(aes(x = n, y = fct_reorder(tax_name, n))) + 
  geom_col() +
  labs(x = 'n distinct RIT elements', y = 'Taxon name',
       title = 'RITs elements by taxon')

distinct_rits |>
  dplyr::count(phylum, class) |> 
  ggplot(aes(x = n, y = fct_reorder(class, n))) + 
  geom_col(show.legend = F) +
  scale_y_discrete(position = 'left') +
  facet_grid(phylum~., scales = 'free_y', space = 'free')   +
  theme(strip.text.y = element_text(angle = 0)) +
  labs(x = 'n distinct RIT elements', y = NULL,
       title = 'RITs elements by phylum, class')

distinct_rits |>
  dplyr::count(phylum, class) |> 
  ggplot(aes(x = n, y = fct_reorder(class, n), fill = phylum)) + 
  geom_col(show.legend = F) +
  scale_y_discrete(position = 'left') +
  facet_grid(phylum~., scales = 'free_y', space = 'free')   +
  theme(strip.text.y = element_text(angle = 0)) +
  labs(x = 'n distinct RIT elements', y = NULL,
       title = 'RITs elements by phylum, class')


## Plot frequency histogram of unique RIT sequences
# 420 unique sequences; some as many as 16 copies within a strain
rit_elements |>
  dplyr::count(tax_id, rit_unique_seq_id) |>
  arrange(desc(n)) |>
  ggplot(aes(n, fill = rit_unique_seq_id)) +
  geom_histogram(position = 'stack', color = 'darkgray', size = 0.15,
                 show.legend = F) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  labs(y = 'n distinct RIT elements', x = 'copy number',
       title = 'Copy number distibution for RIT elements (identical copies)')


# how large are RIT ABC proteins?
rit_elements |> 
  select(protein_df) |> 
  unnest(protein_df) |> 
  filter(str_detect(prot_pos, 'p')) |>
  mutate(prot_length = (prot_stop - prot_start +1)/3) |> 
  ggplot(aes(prot_length, fill = prot_pred)) + 
  geom_histogram() +
  facet_wrap(~prot_pos, ncol = 1, scales = 'free_y') +
  rcartocolor::scale_fill_carto_d() +
  theme_minimal()

# what are the integrase classifier predictions for the flanking genes?
rit_elements |> 
  select(protein_df) |> 
  unnest(protein_df) |> 
  filter(!str_detect(prot_pos, 'p')) |>
  mutate(
    prot_pos = ifelse(prot_pos == 'lflank', 'Left flank', 'Right flank'),
    prot_pred = case_when(
      is.na(prot_pred) ~ 'Pseudogene',
      prot_pred == 'Other' ~ 'Non-integrase',
      TRUE ~ prot_pred),
    prot_pred = fct_rev(fct_relevel(prot_pred, 
                                    c('Non-integrase', 'Pseudogene'), 
                                    after = 0L))
  ) |>
  ggplot(aes(y = prot_pred)) + 
  geom_bar() +
  facet_wrap(~prot_pos, ncol = 1) +
  labs(y = NULL, 
       title = 'What are the genes are adjacent to RIT elements?')+
  theme_minimal()

## How are the adjacent genes annotated?
rit_elements |> 
  select(protein_df) |> 
  unnest(protein_df) |> 
  filter(!str_detect(prot_pos, 'p')) |>
  filter(prot_pred == 'Other') |> 
  mutate(prot_product = 
           str_to_lower(prot_product) |> 
           str_remove(' protein$') |> 
           str_remove('^putative ') |> 
           str_remove(' family$') |> 
           str_remove(' domain-containing$') |> 
           str_remove('-like$') |> 
           str_remove('site-specific dna ') |> 
           str_replace_all('(conserved )?hypothetical.*?|uncharacterised|protein of unknown function', 'unannotated') |> 
           str_replace_all('[Uu]ncharcterised ', 'unannotated') |> 
           str_remove('-containing') |> 
           str_to_lower(),
         prot_product = case_when(
           str_detect(prot_product, 'duf[0-9]+') ~
             str_extract(prot_product, 'duf[0-9]+') |> str_to_upper(),
           str_detect(prot_product, '^dde') ~ 'DDE transposase',
           TRUE ~ prot_product
         )
  ) |> 
  dplyr::count(prot_pos, prot_product) |> 
  mutate(prot_pos = ifelse(prot_pos == 'lflank', 
                           'Left flank', 'Right flank')) |> 
  arrange(desc(n)) |> 
  slice(1:70) |> 
  ggplot(aes(y = fct_reorder(prot_product, n), x = n, fill = prot_pos)) + 
  geom_col() +
  scale_x_log10() +
  rcartocolor::scale_fill_carto_d(palette = 1) +
  labs(y = NULL, fill = NULL, 
       title = 'How are the adjacent genes annotated?') +
  theme_bw()

# 
# 
#   