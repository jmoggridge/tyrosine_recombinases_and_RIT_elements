library(tidyverse)
library(janitor)

aclame <- 
  read_tsv("./Data/aclame_proteins_all_0.4.tab", skip = 1) %>% 
  clean_names()

glimpse(aclame)

aclame_integrase <- aclame %>% 
  filter(str_detect(ncbi_annotation, 'integrase')) %>% 
  filter(!str_detect(ncbi_annotation, 'serine'))

aclame_recombinase <- aclame %>% 
  filter(str_detect(ncbi_annotation, 'tyrosine recombinase'))

aclame_annot_int <- aclame %>% 
  filter(str_detect(aclame_function_s, "tyrosine-based recombinase activity"))
