## RIT distinct results EDA and grouping by similarity
library(tidyverse)

rits <- read_rds('results/non_redund_rits_NJ_clustered.rds') |> 
  mutate(obs_id = row_number())
rits


spirochaetes <- rits |> 
  filter(phylum =='Spirochaetes')

spirochaetes |> View()

taxonomy <- read_rds('data/CDD/tax_data.rds') |> 
  unnest(cols = c(tax_lineage)) |> 
  select(tax_id, tax_name, clade, phylum, class, order, family, genus, species, strain) |> 
  rename_with(~ifelse(str_detect(.x, 'tax_'), .x, paste0('tax_', .x)))


rits |> 
  left_join(taxonomy) |> 
  select(obs_id, tax_id, nuc_name, clade, phylum, class, order, family, genus, species, strain, contains('tax')) |> 
  View()

# 95% similarity clusters: phyla, families, organisms
phyla_clusts <- rits |> 
  group_by(cluster0_05NJ) |> 
  select(cluster0_05NJ, phylum, family, organism) |> 
  summarize(
    phyla = list(unique(phylum)),
    fams = list(unique(family)),
    organisms = list(unique(organism))
    ) |> 
  mutate(n_phyla = map_dbl(phyla, length),
         n_fams = map_dbl(fams, length)) |> 
  arrange(desc(n_phyla)) |> 
  unnest(phyla) |> 
  count(phyla) 
  
old_phyla <- c(
  'Alphaproteobacteria', 'Betaproteobacteria', 'Gammaproteobacteria', 'Delta/Epsilon',
  'Firmicutes', 'Verrucomicrobia', 'Bacteroidetes/Chlorobi', 'Actinobacteria', 'Acidobacteria'
)
phyla_clusts |> 
  rename(phylum  = phyla) |> 
  mutate(phylum = ifelse(phylum == 'Spirochaetia', 'Spirochaetes', phylum)) |> 
  mutate(phylum = ifelse(phylum == 'Unclassified', 'Candidate', phylum)) |> 
  mutate(color = phylum %in% old_phyla) |> 
  left_join(rits |> select(clade, phylum) |> distinct()) |> 
  ggplot(aes(n, fct_rev(phylum), fill = color)) +
  geom_col(show.legend = F) +
  labs(y = 'Phylum', x = 'n 95% similarity clusters') +
  theme_bw() +
  facet_grid(clade~., scales = 'free_y', space = 'free') +
  theme(strip.text.y = element_text(angle = 0))
  
rits |> filter(str_detect(clade, 'Nitrospinae')) |> View()
rits |> 
  group_by(cluster0_05NJ) |> 
  select(cluster0_05NJ, phylum, family)
  # 
  # unnest(phyla) |> 
  # mutate(n = row_number()) |> 
  # ungroup() |> 
  # group_by(phyla) |> 
  # nest(rits = c(rit_id,rit_dna))

  # mutate(n_phyla = map_dbl(phyla, length),
         # n_fam = map_dbl(fams, length),
         # n_org = map_dbl(organisms, length)) |> 
  

distinct_rits <- 
  rits |> 
  group_by(tax_id, rit_id) |> 
  nest(rit_occurrences = c(
    nuc_id, nuc_accession, nuc_name, rit_id,
    rit_length, protein_df, matches('rit_st|upstream|downstream')
    )
  ) |> 
  ungroup() |> 
  mutate(copies = map_int(rit_occurrences, ~nrow(.x))) 

distinct_rits |> View()

# TODO **redo with clustered RITs ***

## EDA taxonomy ============================================================

# 670 distinct taxon ids; one taxon has 35! distinct RITs.
distinct_rits |> dplyr::count(tax_id) |> arrange(desc(n))
distinct_rits |> dplyr::count(tax_id, organism) |> print(n=200)

distinct_rits |>
  dplyr::count(genus, species, strain) |> 
  arrange(desc(n)) |> 
  filter(!is.na(strain)) |> 
  ggplot(aes(n)) + 
  geom_histogram() + 
  labs(
    x = 'number of distinct elements',
    y = 'number of taxa',
    title = 'Distribution of the number of distinct elements per strain'
    )

# one strain has 32 RITs. -> probably just SNPs??
distinct_rits |>
  dplyr::count(phylum, tax_name) |> 
  filter(n > 2) |> 
  ggplot(aes(x = n, 
             y = fct_reorder(tax_name, n), 
             color = phylum,
             fill = phylum)) + 
  geom_col() +
  labs(x = 'n distinct RIT elements', y = 'Taxon name',
       title = 'RITs elements by taxon')

distinct_rits |>
  dplyr::count(phylum, class) |> 
  ggplot(aes(x = n, y = fct_rev(class), fill = phylum)) + 
  geom_col(show.legend = F) +
  scale_y_discrete(position = 'left') +
  facet_grid(phylum~., scales = 'free_y', space = 'free')   +
  theme(strip.text.y = element_text(angle = 0)) +
  labs(x = 'n distinct RIT elements', y = NULL,
       title = 'RITs elements by phylum, class')



## Plot frequency histogram of unique RIT sequences
# 420 unique sequences; some as many as 16 copies within a strain
rits |>
  dplyr::count(tax_id, rit_unique_seq_id) |>
  arrange(desc(n)) |>
  ggplot(aes(n, fill = rit_unique_seq_id)) +
  geom_histogram(position = 'stack', color = 'darkgray', size = 0.15,
                 show.legend = F) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  labs(y = 'n distinct RIT elements', x = 'copy number',
       title = 'Copy number distibution for RIT elements (identical copies)')


# how large are RIT ABC proteins?
rits |> 
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












# P2_AXX Visualizations

rits <- read_rds('results/non_redund_rits_NJ_clustered.rds')

# From 1087 nucleotide records,
# hard to figure out what is partial or complete / contigs v scaffolds
# finished vs draft.
# - ~160 complete genome/plasmid
# - 38 records indicate draft genome
# rest are contig-level sequences.
rits |> 
  select(nuc_id, completeness, sourcedb, genome, nuc_name) |> 
  distinct() |> 
  mutate(
    st = map_chr(
      nuc_name,
      ~{
        tolower(.x) |> 
          str_extract_all(
            'draft|complete genome|wgs|contig|scaffold|complete plasmid'
          ) |> 
          unique() |> 
          paste0(collapse = ' ')
      }),
    st2 = glue('{st}_{completeness}_{genome}_{sourcedb}')
  ) |> 
  # View()
  count(st2) |> 
  print(n=30)


rits |> 
  select(nuc_id, completeness, sourcedb, slen, genome, nuc_name) |> 
  distinct() |> 
  pull(slen) |> 
  summary()

# can see that 'chromosome' mostly genomes but there are far more contigs ~3kb
rits |> 
  select(nuc_id, completeness, sourcedb, slen, genome, nuc_name) |> 
  distinct() |> 
  ggplot(aes(slen, fill = genome)) +
  geom_histogram() +
  scale_x_log10() +
  labs(
    fill = 'ncbi label',
    title = 'how long are the contigs or scaffolds scaffolds?')


rits



## Inverted repeats figures ================================================


rit_IRs_df <- read_rds('results/nr_rits_clustered_IRs.rds')

glimpse(rit_IRs_df)

# length, matches
rit_IRs_df |> 
  pivot_longer(c(al_match, al_percent_id, al_length)) |> 
  ggplot(aes(value)) +
  geom_histogram() +
  facet_wrap(~name, scales = 'free') +
  theme_bw() +
  labs(title = 'Inverted repeat alignments for active RITs (multiple-copy per contig)')

rit_IRs_df |> 
  filter(ir_5p_start > -400) |> 
  filter(ir_3p_start < 400) |> 
  pivot_longer(ir_5p_start:ir_3p_end) |> 
  mutate(
    side = ifelse(str_detect(name, '5p'), 'upstream', 'downstream'),
    startstop = ifelse(str_detect(name, 'start'), 'start', 'stop')) |> 
  
  ggplot(aes(value, fill = startstop, color = startstop)) +
  geom_density(alpha =0.2) +
  facet_wrap(~fct_rev(side), scales = 'free') +
  
  theme_bw() +
  labs(
    x = 'Position relative to RitA start / RitC stop',
    subtitle =
      'Where are the aligned inverted repeats located relative to beginning and end \
of the coding sequence?', 
color = NULL, fill = NULL)

# left TIRs

rit_IRs_df |> 
  select(rit_id, tax_id, nuc_id, al_percent_id, ir_5p_start, ir_5p_end, ir_3p_start, ir_3p_end) |> 
  View()

A <- rit_IRs_df |> 
  mutate(row = row_number()) |> 
  ggplot(aes(y = row, yend = row, x = ir_5p_start, xend = ir_5p_end, color = al_percent_id)) +
  geom_segment(show.legend = F) +
  scale_color_viridis_c(option = 'D', direction = -1) +
  theme_bw() +
  labs(y = 'RIT element', x = 'Upstream position')

B <- rit_IRs_df |> 
  mutate(row = row_number()) |> 
  ggplot(aes(y = row, yend = row, x = ir_3p_start, xend = ir_3p_end, color = al_percent_id)) +
  geom_segment() +
  scale_color_viridis_c(option = 'D', direction = -1) +
  theme_bw() +
  labs(y = NULL, x = 'Downstream position')

library(patchwork)
A + B & plot_layout(guides = 'collect')





