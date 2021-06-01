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

# 
# ## The ICEberg2.0 database has a set of 26,566 proteins associated with integrative conjugative genetic elements (ICEs) as those or predicted to be.
# 
# # Save a fasta file with all the protein sequences from iceberg
# writeLines(readLines("https://db-mml.sjtu.edu.cn/ICEberg2/download/ICE_aa_all.fas"), con = file("./Data/ICEberg_all_ICEs_protein.fasta"))
# 
# # Download fasta of all ICEberg protein sequences for filtering
# ice_fasta <- readDNAStringSet("https://db-mml.sjtu.edu.cn/ICEberg2/download/ICE_aa_all.fas")
# 
# # Create a dataframe to filter integrases from all ICEberg proteins
# ice_df <- 
#   tibble(
#     name = names(ice_fasta),
#     seq = paste(ice_fasta),
#     type = "ICE"
#   ) %>% 
#   # extract annotation info from sequence headers
#   mutate(
#     header = map(name, ~unlist(str_split(.x, "\\|"))),
#     database = map_chr(header, 1),
#     gi = map_chr(header, 2),
#     not_sure_what_this_is = map_chr(header, 3),
#     # seqs from genbank, refseq, emb, dbj
#     source = map_chr(header, 4),
#     accession = map_chr(header, 5),
#     annotation = map_chr(header, 6),
#     # which ICE
#     element = str_extract(annotation, "\\[.+\\]"),
#     element = str_remove_all(element, "\\[|\\]"),
#     # protein annotation
#     protein = str_remove(annotation, " \\[.+\\]")
#   ) %>% 
#   select(-annotation)
# glimpse(ice_df)  
# 
# integrase <- ice_df %>% 
#   filter(str_detect(protein, 'integrase'))
# 
# summary(factor(ice_df$source))
# 
# view(
#   ice_df %>% 
#     group_by(protein) %>% 
#     filter()
# )
# 


