# P2_B3 trying to find answers for nicoles question about the two strains with RITs in ICEberg..

library(tidyverse)
library(Biostrings)
library(tidymodels)
library(furrr)
library(patchwork)
library(seqinr)
library(glue)

# full dataset
iceberg <- read_tsv('data/iceberg/ICE_db.tsv') |> 
  select(-parent_element, -prot_mystery_id)
glimpse(iceberg)

results_with_RITs <- read_csv('./results/iceberg_RITs.csv')
glimpse(results_with_RITs)

## PART 2 Examine RIT elements ------

# Nicole wants to know if the Rits in these strains are identical (similarity)...

# reorganize dataset
rits <- 
  results_with_RITs |> 
  left_join(
    iceberg |> 
      select(-type, -id), by = c("name", "accession", "description", 
                            "prot_name", "prot_database", "prot_accession")) |> 
  mutate(across(contains('pred'), as.character)) |> 
  mutate(
    strain = str_remove(description, '[0-9]+..[0-9]+\\s'),
    # keep consensus of predictions or best with prob > 0.5
    prediction = case_when(
      knn_pred %in% glmnet_pred ~ paste(knn_pred),
      knn_prob > glmnet_prob & knn_prob > 0.5 ~ paste(knn_pred),
      knn_prob < glmnet_prob & glmnet_prob > 0.5 ~ paste(glmnet_pred),
      TRUE ~ 'no consenus')
  ) |> 
  select(-c(description, contains('_pred'), contains('_prob'))) |> 
  ungroup() |> 
  # rename for clarity
  transmute(
    strain, 
    ice_id = id, 
    ice_name = name, 
    ice_types = type,
    ice_start = start_pos, 
    ice_stop = end_pos,  
    dna_acc = accession, 
    dna_db = database, 
    dna_seq, 
    dna_length = map_dbl(dna_seq, nchar),
    prot_acc = prot_accession, prot_name, prot_seq,
    prot_length = map_dbl(prot_seq, nchar),
    prediction
  ) |> 
  distinct()


glimpse(rits)

# there are three ICEs from two strains with putative RIT elements...
rits |> dplyr::count(strain, ice_name)

# look and guess which are rit elements
rits |> select(-contains('seq')) |> View()

# created a table with each Rit element as a row with proteins 1 - 3
elements <- tribble(
  ~ice_name,   ~p1_acc,       ~p2_acc,      ~p3_acc,  ~ arrangement, ~comment,
  'ICE-GI1', 'CAP41371.1', 'CAP41372.1', 'CAP41373.1',  'CBA', '',
  'ICE-GI2', 'CAP41638.1', 'CAP41639.1', 'CAP41640.1',  'ABC', '',
  'ICE-GI3', 'CAP41785.1', 'CAP41786.1', 'CAP41789.1',  'ABC', 'contains transposases',
  'ICE-GI3', 'CAP41836.1', 'CAP41837.1', 'CAP41838.1',  'ABC', '', 
  'ICEMlSym(NZP2037)', 'ANN60838.1', 'ANN60839.1', 'ANN60840.1', 'CBA', '',
  'ICEMlSym(NZP2037)', 'ANN60865.1', 'ANN60866.1', 'ANN60867.1', 'ABC', ''
)

glimpse(elements)

# prep the protein sequences to join to elements' p1-3
proteins <- 
  rits |> 
  ungroup() |> 
  filter(prot_acc %in% c(elements$p1_acc, elements$p2_acc, 
                               elements$p3_acc))
p1 <- proteins |> 
  transmute(
    p1_acc = prot_acc,
    p1_seq = prot_seq,
    p1_length = prot_length
    )
p2 <- proteins |> 
  transmute(
    p2_acc = prot_acc,
    p2_seq = prot_seq,
    p2_length = prot_length
  )
p3 <- proteins |> 
  transmute(
    p3_acc = prot_acc,
    p3_seq = prot_seq,
    p3_length = prot_length
  )

# join the protein data to elements, join dna data to elements.
rit_elements <- elements |> 
  left_join(p1) |> left_join(p2) |> left_join(p3) |> 
  distinct()  |> 
  left_join(
    rits |> select(contains('ice'), contains('dna')),
    by = 'ice_name') |> 
  distinct() |> 
  ## fix the name of the ICE from M. loti
  mutate(
    ice_name = ifelse(str_detect(ice_name, 'ICEMlSym'),
                      'ICEMlSym(NZP2037)_alpha',
                      ice_name)
  )

rit_elements |> select(-contains('seq')) |> View()

# rm(proteins, p1, p2, p3, elements, iceberg)

glimpse(rit_elements)


## TRANLATION CODE THAT I GAVE UP ON

# rit_elements |> pull(dna_acc) |> unique()
# rit_elements |> pull(dna_seq) |> unique()
# rit_elements |> pull(dna_length) |> unique()
# 
# rit_elements |> select(-contains('seq')) |>  View()
# 
# # translation codes are all 11, add strains
# rit_elements <- rit_elements |> 
#   left_join(
#     tribble(
#       ~dna_acc, ~dna_code,
#       'AM902716', 11, 
#       'CP016079', 11
#     )) |> 
#   distinct()  |> 
#   left_join(rits |> select(strain, ice_id) |> distinct())
# 
# # translate dna string for all reading frames, code = 11
# make_translations <- function(dna){
#   dna <- s2c(dna)
#   tibble(
#     frame = glue('transl_{1:6}'),
#     offset = rep(c(0,1,2), 2), 
#     orientation = rep(c('F','R'), each = 3),
#     translation = map2_chr(
#       .x = offset,
#       .y = orientation,
#       .f = ~{
#         seqinr::translate(dna, frame = .x, sens = .y, numcode = 11) |> 
#           c2s()
#       }) 
#   )
# }
# 
# # locate proteins in translations.
# locate_protein <- function(translations, protein, tag){
#   find_pos <- 
#     translations |> 
#     mutate(
#       prot_acc = tag,
#       start = str_locate(translation, protein)[,1],
#       end = str_locate(translation, protein)[,2],
#       length = end - start + 1,
#     )
#   return(find_pos)
# }
# 
# 
# #  map proteins to gene locations
# translations <- map(rit_elements$dna_seq, make_translations)
# locations <- pmap(
#   .l = list(translations,  rit_elements$p1_seq, rit_elements$p1_acc), 
#   .f = function(x,y,z) locate_protein(x,y,z)
# )
# locations |> bind_rows() |> select(-translation) |> View()
# 








## NOTE: Couldnt manage to find proteins with ICE translations... records don't seem to match up properly....


# 
# rit_elements |> select(-contains('seq')) |> View()
# 
# rit_elements
# 
# 
# dna <- rit_elements$dna_seq[1]
# translations <- make_translations(dna)
# 
# dna_length <- nchar(dna)
# prot_length <- nchar(protein)
# dna_start <- (location$start-1)*3 + 1 + location$offset
# dna_end <- (location$end-1)*3 + 3 + location$offset
# 
# 
# ((dna_end - dna_start + 1)/3) == prot_length
# 
# 
# 
# dna_start_actual <- ifelse(
#   location$orientation == 'R',
#   yes = dna_length - dna_end+1,
#   no = dna_start
#   )
# 
# dna_end_actual <- ifelse(
#   location$orientation == 'R',
#   yes = dna_length - dna_start+1,
#   no = dna_end
# )
# 
# dna_c <- s2c(dna)
# dna_segment <- ifelse(
#   location$orientation == 'R',
#   yes = c2s(rev(comp(dna_c[dna_start_actual: dna_end_actual]))),
#   no = c2s(dna_c[dna_start_actual: dna_end_actual] )
# )
# dna_segment
# c2s(translate(s2c(dna_segment))) == protein
# 
