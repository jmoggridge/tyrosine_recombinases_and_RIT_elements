library(tidyverse)

genomes <- read_csv('data/genome_info_by_organism.csv') |> 
  janitor::clean_names() |> 
  separate(organism_groups, ';', into = c('domain', 'clade', 'phylum')) |> 
  mutate(across(where(is.character), as_factor)) |> 
  mutate(
    phylum = str_remove(phylum, ' group| subdivisions') |> 
      str_trim() |> 
      str_to_title(),
    phylum = ifelse(phylum == 'Acidobacteriia', 'Acidobacteria', phylum),
  ) 
  # filter(!clade %in% c('unclassified Bacteria', 'Bacteria incertae sedis')) |> 
  # filter(!str_detect(phylum, '^Candidat|^Enviro'))

genome_size <- genomes |> 
  group_by(phylum) |> 
  summarize(median_size = median(size_mb),
            max_size = max(size_mb),
            min_size = min(size_mb),
            q1 = quantile(size_mb, 0.25),
            q3 = quantile(size_mb, 0.75)
            )
# need to add a row for stenosarchaea


rits <- read_rds('results/non_redund_rits_NJ_clustered.rds') |> 
  select(phylum) |> 
  mutate(
    phylum = ifelse(phylum == 'Spirochaetes', 'Spirochaetia', phylum),
    phylum = ifelse(phylum == 'Acidobacteriia', 'Acidobacteria', phylum),
    phylum = str_to_title(phylum) |> str_trim(),
    obs_id = row_number()
  ) |> 
  filter(!str_detect(phylum, '^Candidat|^Enviro')) |> 
  mutate(across(where(is.character), as_factor)) 
  

genomes |> summary()
rits |> summary()

df <- 
  rits |> 
  dplyr::count(phylum, name = 'RITs') |>
  arrange(desc(RITs)) |> 
  full_join(
    genomes |> dplyr::count(phylum, name = 'NCBI') |> arrange(desc(NCBI))
  ) |>
  full_join(genome_size) |> 
  arrange(desc(RITs)) |> 
  filter(RITs > 0 & !phylum %in% c('Unclassified', 'Candidate', 'Stenosarchaea')) |>
  mutate(across(everything(), ~replace_na(.x, 0))) |> 
  mutate(NCBI = NCBI / sum(NCBI),
         RITs = RITs / sum(RITs))

A <- df |> 
  pivot_longer(cols = c(NCBI, RITs), names_to = 'group', values_to = 'prop') |> 
  ggplot(aes(y = fct_rev(phylum), x = prop, fill = group)) +
  geom_col(position = 'dodge', width = 0.8, alpha = 0.8) +
  labs(x = 'proportion of group', fill = 'Group', y = '') +
  theme(axis.text.y.left = element_text(size = 5)) +
  theme_bw()


B <- df |> 
  ggplot(aes(y = fct_rev(phylum))) +
  geom_pointrange(aes(x = median_size, xmin = q1, xmax = q3), fatten = 2) +
  geom_point(aes(x = max_size), shape = 1, show.legend = T) +
  geom_point(aes(x = min_size), shape = 1, show.legend = T) +
  scale_x_log10() +
  labs(y = NULL, x = 'Genome size (Mb)') +
  theme_bw() +
  theme(axis.text.y.left = element_blank(),
        axis.ticks.y = element_blank())
  

library(patchwork)
A + B + plot_layout(guides = 'collect')

