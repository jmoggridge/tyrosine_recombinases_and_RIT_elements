## RIT Data for Nicole

library(tidyverse)

rit_elements <- read_rds('./results/rit_elements.rds') |> 
  left_join(read_rds('./data/CDD/tax_data_fixed.rds')) |> 
  unnest(tax_lineage)

# 185 unique taxa
rit_elements |> 
  pull(tax_id) |> 
  unique() |> 
  length()

firmi_gamma <- rit_elements |> 
  filter(phylum %in% c('Gammaproteobacteria', 'Firmicutes')) |> 
  mutate(nuc_name = map2_chr(nuc_name, nuc_accession, 
                         ~str_remove(.x, .y) |> str_trim())) |> 
  select(tax_id, tax_name, nuc_accession, 
         rit_id, rit_start, rit_stop, rit_strand) |> 
  mutate(tax_id = as.numeric(tax_id)) |> 
  arrange(tax_id, tax_name, nuc_accession)

firmi_gamma |> print(n=150)
write_csv(firmi_gamma, 'RITs_firmicutes_gammaproteobacteria.csv')

  