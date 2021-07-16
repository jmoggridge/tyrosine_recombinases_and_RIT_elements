## P2_A5 figure out how to identify rit element from genbank...

library(tidyverse)
library(progress)
library(glue)
library(rentrez)
source('./R/P2_parse_genbank_xml.R')
set_entrez_key('889fdb9786a14019a0a1257196a09ba4ba08')
Sys.getenv('ENTREZ_KEY')

## These genbank records have no CDS included.... problem!
# no_CDS <- c('1834417989', '1024364864')

three_ints <- read_rds('./data/CDD/ids_w_three_integrases.rds') |> 
  arrange(slen)
three_ints


# for one nucleotide -> get genbank -> identify RITs



# 
# 
# ## PARSE FILES -------
# rm(count_table, id_data, nuc_summary, three_integrases, to_retrieve)
# 
# wont_parse <- c()
# 
# gb1 <- read_rds('./data/CDD/RIT_gbk_1.rds') |> slice(1:20)
# 
# pb <- progress_bar$new(total = 20)
# parsed_gb1a <- gb1 |> 
#   filter(!nuc_id %in% wont_parse) |> 
#   filter(!nuc_id %in% no_CDS) |> 
#   mutate(genbank_rcd = map2(gbk, nuc_id, ~{
#     print(glue('\nnuc_id: {.y}'))
#     parsed <- parse_genbank(.x)
#     pb$tick()
#     return(parsed)
#   })
#   ) |> 
#   select(nuc_id, genbank_rcd)
# write_rds(parsed_gb1a, './data/CDD/RIT_gbk_1a_parsed.rds')
# beepr::beep()
# rm(parsed_gb1a)
# 
# 
# gb1 <- read_rds('./data/CDD/RIT_gbk_1.rds') |> slice(21:40)
# pb <- progress_bar$new(total = 20)
# parsed_gb <- gb1 |> 
#   filter(!nuc_id %in% wont_parse) |> 
#   filter(!nuc_id %in% no_CDS) |> 
#   mutate(genbank_rcd = map2(gbk, nuc_id, ~{
#     print(glue('\nnuc_id: {.y}'))
#     parsed <- parse_genbank(.x)
#     pb$tick()
#     return(parsed)
#   })
#   ) |> 
#   select(nuc_id, genbank_rcd)
# write_rds(parsed_gb_1b, './data/CDD/RIT_gbk_1b_parsed.rds')
# beepr::beep()
# 
# 
# 
# 
# 
# 
