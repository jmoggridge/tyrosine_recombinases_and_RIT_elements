## P2_A5 figure out how to identify rit element from genbank...

library(tidyverse)
library(progress)
library(glue)
library(rentrez)
source('./R/P2_parse_genbank_xml.R')
set_entrez_key('889fdb9786a14019a0a1257196a09ba4ba08')
Sys.getenv('ENTREZ_KEY')
# genbank files for these nuc ids lack the CDS features.
no_cds <- c('1834417989', '1024364864')

# these nucleotides probably contain RIT element bc have 3 CDD proteins.
three_ints <- read_rds('./data/CDD/ids_w_three_integrases.rds') |> 
  filter(!nuc_id %in% no_cds) |> 
  arrange(slen)
three_ints

# prot_id and accessions for this set of nt's.
three_ints_prots <- three_ints |> unnest(prot_id) |> pull(prot_id)

# join with other accession numbers
prot_accessions <- read_rds('./data/CDD/prot_summary_fixed.rds') |> 
  filter(prot_id %in% three_ints_prots) |> 
  select(prot_id, caption)

# create an index for nuc ids to their genbank records
genbank_files_index <- 
  tibble(file = Sys.glob('./data/CDD/RIT_gbk_[0-9].rds')) |> 
  mutate(nuc_id = map(file, ~ read_rds(.x) |> pull(nuc_id))) |> 
  unnest(nuc_id)



# for a nuc_id, return the parsed genbank record as a tibble
open_genbank <- function(x){
  # figure out which file to open from index
  file <- genbank_files_index |> 
    filter(nuc_id == x) |> 
    pull(file)
  file
  
  # pull the genbank record and parse
  gbk <- read_rds(file) |> 
    filter(nuc_id == x) |> 
    unnest(gbk) |> 
    pull(gbk) |> 
    parse_genbank()
  return(gbk)
}

# example: first nuc_id
x <- three_ints$nuc_id[[1]]
gbk <- open_genbank(x)

# are the three RIT proteins there??
prot_x <- three_ints |>
  filter(nuc_id == x) |> 
  unnest(prot_id) |> pull(prot_id)
cap_x <- prot_summary |> 
  filter(uid %in% prot_x) |> 
  pull(caption) |>
  map_chr(~glue('{.x}.1'))

# extract features table
ft <- gbk |>
  select(feature_table) |> 
  unnest(feature_table) |> 
  select(-gb_feature_key, -locus_tag)
ft |> 
  View()


which(ft$protein_id %in% cap_x)


# classify proteins

classify_proteins <- function(gbk){
  
  
}

## HMM scores
# make temp fasta file of seqs to score against hmm library
fasta_path <- tempfile()
fasta <- Biostrings::AAStringSet(gbk$translation)
names(fasta) <- gbk${{protein}}
writeXStringSet(fasta, fasta_path)  



# figure out where there is three RITs in a row...


gbk_parsed <- parse_genbank(gbk)
gbk_parsed
feature_table <- gbk_parsed |>
  select(feature_table) |> 
  unnest(feature_table)
gbk_prot_id <- feature_table |> pull(protein_id)

prot_for_nuc_id <- 
  three_ints |> 
  filter(nuc_id == x) |> 
  pull(prot_id)
any(gbk_prot_id)


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
