library(tidyverse)
library(janitor)
library(glue)
library(scales)

# genomes 
genome_list <- 
  read_csv("./Data/ncbi_prokaryotes_genomes_list.csv") %>% 
  clean_names() %>% 
  mutate(type = "Genome")
glimpse(genome_list)

# plasmids
plasmid_list <- 
  read_csv("./Data/ncbi_plasmids_list.csv") %>% 
  clean_names() %>% 
  mutate(type = "Plasmid")
glimpse(plasmid_list)

# combine data and plot cumulative distribution over time
bind_rows(
  genome_list %>% select(release_date, type),
  plasmid_list %>% select(release_date, type)
) %>% 
  ggplot(aes(release_date, color = type)) +
  stat_ecdf(pad = F) +
  theme_light() +
  theme(legend.position = c(0.1, 0.85)) +
  labs(x = NULL, y = 'Proportion', color = NULL,
       subtitle = glue("Cumulative distribution of release dates of {comma(nrow(genome_list))} bacterial genomes \nand {comma(nrow(plasmid_list))} plasmids in NCBI Genome database (2021-05-19)"))

