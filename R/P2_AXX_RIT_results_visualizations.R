## RIT distinct results EDA and grouping by similarity
library(tidyverse)
library(patchwork)
library(glue)

## Data ---------------------------------------------------------------------
rits <- read_rds('results/non_redund_rits_NJ_clustered.rds') |> 
  mutate(
    phylum = ifelse(phylum == 'Spirochaetia', 'Spirochaetes', phylum),
    obs_id = row_number()
    ) 

rits |> 
  unnest(protein_df) |> 
  select(-contains('dna'), -contains('seq')) |> 
  filter(prot_pos %in% c('lflank', 'rflank')) |> 
  View()


rits |> 
  count(rit_id, nuc_id) |> 
  arrange(desc(n))
unique(rits$genus)

rits
length(unique(rits$nuc_id))
length(unique(rits$phylum))

# Are all RITs RitABC arrangement?
rits |> 
  unnest(protein_df) |> 
  select(prot_pos, prot_pred) |> 
  filter(prot_pos %in% c('p1', 'p2', 'p3')) |> 
  pivot_wider(names_from = prot_pos,
              values_from = prot_pred, 
              values_fn = list) |> 
  unnest() |> 
  mutate(across(everything(), as_factor)) |> 
  summary()

# RIT length histogram
rits |> 
  select(rit_length) |> 
  ggplot() +
  geom_histogram(aes(x = rit_length)) +
  theme_bw() +
  labs(x = 'RIT coding sequence length')

# Some of these really short ones might be partial elements??
rits |> 
  filter(rit_length < 2500) |> 
  View()

# active_rits <- read_rds('results/')

# Gene overlaps -------------------------------------------------------------
rits |> 
  select(rit_id, nuc_id, protein_df) |> 
  unnest(protein_df) |> 
  relocate(prot_overlap) |> 
  filter(prot_pos %in% c('p1', 'p2')) |> 
  select(prot_pos, prot_overlap) |> 
  pivot_wider(names_from = prot_pos, 
              values_from = prot_overlap, 
              values_fn = list) |> 
  unnest() |> 
  summary()

# some very large gaps / overlaps
# rits |> 
#   select(rit_id, nuc_id, protein_df) |> 
#   unnest(protein_df) |> 
#   relocate(prot_overlap) |> 
#   filter(prot_pos %in% c('p1', 'p2')) |> 
#   filter(prot_overlap < -100) |> 
#   View()

# overlaps histograms
rits |> 
  select(rit_id, nuc_id, protein_df) |> 
  unnest(protein_df) |> 
  relocate(prot_overlap) |> 
  filter(prot_pos %in% c('p1', 'p2')) |> 
  mutate(prot_pos = ifelse(prot_pos == 'p1', 'RitA-RitB', 'RitB-RitC')) |> 
  select(prot_pos, prot_overlap) |> 
  ggplot(aes(prot_overlap)) +
  geom_histogram() +
  facet_wrap(~prot_pos) +
  labs(x = 'Gap (-) / Overlap (+)') +
  theme_bw()

## Spirochaetes ----------------------------------------------------------------
# Dr. Kropinski was interested in the Spirochaetes RIT...
spirochaetes <- rits |> 
  filter(tolower(phylum) %in% c('spirochaetes', 'spirochaetia'))
spirochaetes |> View()

## The taxonomy labels need fixing in many cases - 
# strain and species frequently missing their labels

# taxonomy <- read_rds('data/CDD/tax_data.rds') |> 
#   unnest(cols = c(tax_lineage)) |> 
#   select(tax_id, tax_name, clade, phylum, class, order, family, genus, species, strain) |> 
#   rename_with(~ifelse(str_detect(.x, 'tax_'), .x, paste0('tax_', .x)))
# 
# # check taxonomy labels
# rits |> 
#   left_join(taxonomy) |> 
#   select(obs_id, tax_id, nuc_name, clade, phylum, class, order, family, genus, species, strain, contains('tax')) |> 
#   View()

## Visualize taxonomy 1 --------------------------------------------------------

# phyla with prev. described RITs
old_phyla <- c(
  'Alphaproteobacteria', 'Betaproteobacteria', 'Gammaproteobacteria',
  'Delta/Epsilon','Firmicutes', 'Verrucomicrobia', 'Bacteroidetes/Chlorobi',
  'Actinobacteria', 'Acidobacteria')

# number of 95% similarity clusters: phyla, families, organisms
phyla_clust_95 <- rits |> 
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
  
phyla_clust_95 |> 
  mutate(phylum = phyla) |> 
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
  

## For 70% NJ clusters:
phyla_clust_70 <- rits |> 
  group_by(cluster0_3NJ) |> 
  select(cluster0_3NJ, phylum, class, organism) |> 
  summarize(
    phyla = list(unique(phylum)),
    class = list(unique(class)),
    organisms = list(unique(organism))
  ) |> 
  mutate(n_phyla = map_dbl(phyla, length),
         n_class = map_dbl(class, length)) |> 
  arrange(desc(n_class)) |> 
  unnest(class) |> 
  count(class)  |> 
  left_join(rits |> select(phylum, class) |> distinct()) |> 
  mutate(phylum = ifelse(phylum == 'Spirochaetia', 'Spirochaetes', phylum),
         phylum = ifelse(phylum == 'Unclassified', 'Candidate', phylum),
         color = phylum %in% old_phyla,
         class = str_replace(class, 'unclassified Unclassified', 'Unclassified'),
         class = str_to_title(class),
         ) |> 
  arrange(phylum)
  

phyla_clust_70 |> 
  ggplot(aes(n, fct_rev(class), fill = color)) +
  geom_col(show.legend = F) +
  labs(y = NULL, x = 'number of 70% similarity clusters (by NJ)') +
  theme_bw() +
  facet_grid(phylum~., scales = 'free_y', space = 'free') +
  theme(strip.text.y = element_text(angle = 0),
        panel.spacing = unit(0, 'mm')
        )

##
## 70% NJ clusters
phyla_clust_50 <- rits |> 
  group_by(cluster0_5NJ) |> 
  select(cluster0_5NJ, phylum, class, organism) |> 
  summarize(
    phyla = list(unique(phylum)),
    class = list(unique(class)),
    organisms = list(unique(organism))
  ) |> 
  mutate(n_phyla = map_dbl(phyla, length),
         n_class = map_dbl(class, length)) |> 
  arrange(desc(n_class)) |> 
  unnest(class) |> 
  count(class)  |> 
  left_join(rits |> select(phylum, class) |> distinct()) |> 
  mutate(phylum = ifelse(phylum == 'Spirochaetia', 'Spirochaetes', phylum)) |> 
  mutate(phylum = ifelse(phylum == 'Unclassified', 'Candidate', phylum)) |> 
  mutate(color = phylum %in% old_phyla) |> 
  arrange(phylum)

phyla_clust_50 |> 
  # slice(1:91) |> 
  ggplot(aes(n, fct_rev(class), fill = color)) +
  geom_col(show.legend = F) +
  labs(y = 'Phylum', x = 'n 50% similarity clusters (by NJ)') +
  theme_bw() +
  facet_grid(phylum~., scales = 'free_y', space = 'free') +
  theme(strip.text.y = element_text(angle = 0),
        panel.spacing = unit(0, 'mm')
  )


### ------------------------------------------------------------------------

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

# distinct_rits |> View()

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
  dplyr::count(tax_id, rit_id) |>
  arrange(desc(n)) |>
  ggplot(aes(n, fill = rit_id)) +
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
  mutate(prot_length = (prot_stop - prot_start +1)/3,
         prot_pos = case_when(
           prot_pos == 'p1' ~ 'RitA',
           prot_pos == 'p2' ~ 'RitB',
           prot_pos == 'p3' ~ 'RitC',
         )) |> 
  ggplot(aes(prot_length, fill = prot_pred)) + 
  geom_histogram(show.legend = F) +
  facet_wrap(~prot_pos, ncol = 1, scales = 'free_y') +
  rcartocolor::scale_fill_carto_d() +
  labs(x = 'Protein sequence length') +
  theme_bw()

# what are the integrase classifier predictions for the flanking genes?
A <- rits |> 
  select(protein_df) |> 
  unnest(protein_df) |> 
  filter(!str_detect(prot_pos, 'p')) |>
  mutate(
    prot_pos = ifelse(prot_pos == 'lflank', 'Left flank', 'Right flank'),
    prot_pred = case_when(
      is.na(prot_pred) ~ 'Pseudogene or seq. end',
      prot_pred == 'Other' ~ 'Non-integrase',
      prot_pred == 'no consensus' ~ 'No Consensus',
      TRUE ~ prot_pred),
    prot_pred = fct_rev(fct_relevel(prot_pred, 
                                    c('Non-integrase',
                                      'Pseudogene or seq. end',
                                      'No Consensus'), 
                                    after = 0L))
  ) |>
  ggplot(aes(y = prot_pred)) + 
  geom_bar() +
  facet_wrap(~prot_pos, ncol = 1) +
  labs(y = NULL, 
       # title = 'What are the genes are adjacent to RIT elements?'
       )+
  theme_minimal()

## How are the adjacent genes annotated?
B <- rits |> 
  select(protein_df) |> 
  unnest(protein_df) |> 
  filter(!is.na(prot_product)) |> 
  filter(!str_detect(prot_pos, 'p')) |>
  filter(prot_pred == 'Other') |>
  mutate(prot_product = 
           str_to_lower(prot_product) |> 
           str_remove(' protein$') |> 
           str_remove(' c-terminal') |> 
           str_remove('^putative ') |> 
           str_remove(' family$') |> 
           str_remove(' family') |> 
           str_remove(' domain-containing$') |> 
           str_remove('-like$') |> 
           str_remove('site-specific dna ') |> 
           str_replace_all('(conserved )?hypothetical.*?|uncharacterised|protein of unknown function', 'unannotated') |> 
           str_replace_all('[Uu]ncharcterised ', 'unannotated') |> 
           str_remove('-containing') |> 
           str_to_lower() |> 
           str_replace('transposase is66', 'is66 transposase'),
         prot_product = case_when(
           str_detect(prot_product, 'duf[0-9]+') ~
             str_extract(prot_product, 'duf[0-9]+') |> str_to_upper(),
           str_detect(prot_product, '^dde') ~ 'DDE transposase',
           TRUE ~ prot_product
         )
  ) |> 
  group_by(prot_product) |> 
  summarize(lflank = sum(prot_pos == 'lflank'),
            rflank = sum(prot_pos == 'rflank')) |> 
  mutate(n = lflank +rflank) |> 
  arrange(desc(n)) |> 
  dplyr::slice(1:20) |> 
  ggplot(aes(y = fct_reorder(prot_product, n), x = n)) + 
  geom_col() +
  scale_x_log10() +
  labs(y = NULL, fill = NULL, x = 'count',
       # subtitle = 'Non-integrase genes'
       ) +
  theme_bw()

A + B

## What's this strawberry notch protein??
# rits |> 
#   select(protein_df) |> 
#   unnest(protein_df) |> 
#   filter(!str_detect(prot_pos, 'p')) |>
#   filter(prot_pred == 'Other') |> 
#   filter(str_detect(tolower(prot_product), 'strawberry notch')) |> 
#   View()
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
  mutate(genome = ifelse(is.na(genome), 'no information', genome),
         genome = fct_relevel(genome, 'no information', after = 3)) |> 
  
  ggplot(aes(slen, fill = genome)) +
  geom_histogram() +
  scale_x_log10() +
  rcartocolor::scale_fill_carto_d(palette = 2) +
  labs(
    x = 'DNA sequence length (bp)',
    fill = 'NCBI molecule-type label',
    title = 'How long are the contigs/scaffolds?') +
  theme_bw() +
  theme(legend.position = c(0.75, 0.75))


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
  labs(y = 'RIT element id', x = 'Upstream position', fill = "% identity")

B <- rit_IRs_df |> 
  mutate(row = row_number()) |> 
  ggplot(aes(y = row, yend = row, x = ir_3p_start, xend = ir_3p_end, color = al_percent_id)) +
  geom_segment() +
  scale_color_viridis_c(option = 'D', direction = -1) +
  theme_bw() +
  labs(y = NULL, x = 'Downstream position', color = "% identity")

library(patchwork)
A + B & plot_layout(guides = 'collect')




## Active RITs


## better IRs analysis  ----------------------------------------------------

# rits with alignment, keep 1 representative from each of the 5% similarity clusters
rit_IRs_clustreps <- 
  read_rds('results/nr_rits_clustered_IRs.rds') |> 
  group_by(cluster0_05NJ) |> 
  arrange(desc(al_score)) |> 
  dplyr::slice(1) |> 
  ungroup()

rit_IRs_clustreps |> 
  select(starts_with('al_')) |> 
  summary()


rit_IRs_clustreps |> 
  filter(al_length > 200) |> View()


# length, matches
rit_IRs_clustreps |> 
  pivot_longer(c(al_length, al_percent_id,), names_to = 'name') |> 
  mutate(name = case_when(
    # str_detect(name,'al_match') ~ 'Matching positions',
    name == 'al_length' ~ 'Alignment length',
    name == 'al_percent_id' ~ '% Identity',
    # name == 'al_score' ~ 'Alignment score',
    TRUE ~ 'NA'
  )) |> 
  ggplot(aes(value)) +
  geom_histogram(bins = 40) +
  facet_wrap(~name, scales = 'free', nrow = 1) +
  theme_bw() +
  labs(
    x = NULL,
    title = 'Inverted repeat alignments for active RITs (multiple-copy per contig)')


rit_IRs_clustreps |> 
  pivot_longer(ir_5p_start:ir_3p_end) |> 
  mutate(
    side = ifelse(str_detect(name, '5p'), 'upstream', 'downstream'),
    startstop = ifelse(str_detect(name, 'start'), 'start', 'stop')) |> 
  
  ggplot(aes(value, fill = startstop, color = startstop)) +
  geom_density(alpha =0.2) +
  facet_wrap(~fct_rev(side), scales = 'free') +
  theme_light() +
  labs(
    x = 'bp',
    subtitle =
      'Where are the aligned inverted repeats located relative to beginning and end \
of the coding sequence?', 
color = NULL, fill = NULL)

# left TIRs

# rit_IRs_df |> 
#   select(rit_id, tax_id, nuc_id, al_percent_id, ir_5p_start, ir_5p_end, ir_3p_start, ir_3p_end) |> 
#   View()

A <- rit_IRs_clustreps |> 
  mutate(row = row_number()) |> 
  ggplot(aes(y = row, yend = row, x = ir_5p_start, xend = ir_5p_end, color = al_percent_id)) +
  geom_segment(show.legend = F) +
  scale_color_viridis_c(option = 'D', direction = -1) +
  theme_bw() +
  labs(y = 'RIT element index', x = 'Upstream position')

B <- rit_IRs_clustreps |> 
  mutate(row = row_number()) |> 
  ggplot(aes(y = row, yend = row, x = ir_3p_start, xend = ir_3p_end, color = al_percent_id)) +
  geom_segment() +
  scale_color_viridis_c(option = 'D', direction = -1) +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  labs(y = NULL, x = 'Downstream position', color = "% Complementarity")

A + B & plot_layout(guides = 'collect')


# weird relationship between aligned IR length and % identity
rit_IRs_clustreps |> 
  ggplot(aes(al_length, al_percent_id, color = al_score)) +
  geom_point() +
  scale_color_viridis_c() +
  theme_bw()


# show some sample alignments....

rit_IRs_clustreps$al_align[[1]]
rit_IRs_clustreps$al_align[[25]]
rit_IRs_clustreps$al_align[[46]]

rit_IRs_clustreps |> 
  filter(al_score == max(al_score)) |> 
  relocate(starts_with('al_')) |> 
  View()




