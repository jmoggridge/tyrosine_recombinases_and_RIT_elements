## Predict integrases and find RITs in the ICEberg2.0 dataset ##


library(tidyverse)
library(Biostrings)
library(tidymodels)
library(furrr)
library(patchwork)
library(seqinr)

source('./R/00_functions.R')

dir.create('./data/iceberg/hmm_scores')
## Main --------------------------------------------------------------------

# classifiers & hmms
knn_model <- read_rds('./models/knn_classifier.rds')
glmnet_model <- read_rds('./models/glmnet_classifier.rds')
rf_model <- read_rds('./models/rf_classifier.rds')
hmm_folder <- './models/hmm/'

# full dataset
iceberg <- read_tsv('data/iceberg/ICE_db.tsv') |> 
  select(-parent_element, -prot_mystery_id)
glimpse(iceberg)

# ice berg proteins for classifying
ice_proteins <- iceberg |> 
  transmute(acc = prot_accession, prot_seq) |> 
  filter(!is.na(prot_seq)) |> 
  distinct()
ice_proteins


# HMM scores

# make temp fasta file of seqs to score against hmm library
fasta_path <- tempfile()
fasta <- Biostrings::AAStringSet(ice_proteins$prot_seq)
names(fasta) <- ice_proteins$acc
writeXStringSet(fasta, fasta_path)  


# setup paths for hmms and table outputs; 
# build hmmsearch calls (but not executed until next step)
temp_dir <- tempdir()
junk <- tempfile()
hmmsearches <- 
  tibble(hmm_path = Sys.glob(glue('{hmm_folder}*'))) |> 
  mutate(
    out_path = hmm_path |> 
      str_replace(glue('{hmm_folder}'), glue('./data/iceberg/hmm_scores/')) |> 
      str_replace('\\.hmm', '.tbl'),
    calls = glue('hmmsearch --noali -o {junk} --tblout {out_path} {hmm_path} {fasta_path}')
    )

hmmsearches

# do hmmsearch calls in parallel
plan(multisession, workers = availableCores())
hmmsearches <- hmmsearches |> 
  mutate(hmmsearches = future_map_dbl(calls, ~ system(.x))) 
print(all(hmmsearches$hmmsearches == 0))

join_hmmsearches2 <- function(df, files){
  # read and combine hmmsearch data
  searches <- map(files, read_hmmsearch) |> 
    purrr::reduce(.f = full_join, by = 'acc') 
  # join scores to input data & replace any NAs with zeros
  df |> 
    left_join(searches, by = 'acc') |>
    mutate(across(Arch1:Xer, ~replace_na(.x, 0))) 
}

# get scored data
ice_scored <- 
  join_hmmsearches2(df = ice_proteins,
                    files = hmmsearches$out_path)
ice_scored

# make predictions with 3 models
knn_preds <- predict(knn_model, new_data = ice_scored) |> 
  transmute(knn = .pred_class)
glmnet_preds <- predict(glmnet_model, new_data = ice_scored) |> 
  transmute(glmnet = .pred_class)
rf_preds <- predict(rf_model, new_data = ice_scored) |> 
  transmute(rf = .pred_class)

# get prediction probablities with 3 models
knn_prob <- predict(
  knn_model, new_data = ice_scored, type = 'prob'
  ) |>
  bind_cols(ice_scored |> select(acc)) |>
  pivot_longer(.pred_Arch1:.pred_Other) |>
  group_by(acc) |>
  summarize(knn_prob = max(value))

glmnet_prob <- predict(
  glmnet_model, new_data = ice_scored, type = 'prob'
  ) |>
  bind_cols(ice_scored |> select(acc)) |>
  pivot_longer(.pred_Arch1:.pred_Other) |>
  group_by(acc) |>
  summarize(glmnet_prob = max(value))

rf_prob <- predict(
  rf_model, new_data = ice_scored, type = 'prob'
  ) |>
  bind_cols(ice_scored |> select(acc)) |>
  pivot_longer(.pred_Arch1:.pred_Other) |>
  group_by(acc) |>
  summarize(rf_prob = max(value))


# join predictions
iceberg_preds <- 
  ice_scored |> 
  bind_cols(knn_preds, glmnet_preds, rf_preds) |> 
  left_join(knn_prob) |>
  left_join(glmnet_prob) |> 
  left_join(rf_prob) |> 
  left_join(ice_proteins) |> 
  mutate(prot_accession = acc) |> 
  select(-acc) |> 
  relocate(prot_accession, contains('knn'), contains('glmnet'),
           contains('rf'))
  
# add classifications back to full dataset
iceberg_w_preds <- iceberg |> 
  left_join(iceberg_preds, by = c("prot_seq", "prot_accession"))

# see how predictions match
iceberg_preds |>  dplyr::count(knn, glmnet, rf)
 
# y-axis breaks for barcharts 
ynames <- unique(as.character(iceberg_w_preds$knn))

# 
a <- iceberg_w_preds |> 
  dplyr::count(knn) |> 
  ggplot(aes(x = n, y = fct_rev(knn))) +
  geom_col() +
  geom_text(aes(label = n), nudge_x = 800) +
  theme_bw() +
  labs(x = 'count', y = 'YR subfamily prediction') +
  lims(x = c(0, 18000)) +
  scale_y_discrete(breaks = ynames, labels = ynames, drop = T)
a
# 
b <- iceberg_w_preds |> 
  dplyr::count(glmnet) |> 
  ggplot(aes(x = n, y = fct_rev(glmnet))) +
  geom_col() +
  geom_text(aes(label = n), nudge_x = 800) +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank()) +
  labs(x = 'count') +
  lims(x = c(0, 18000)) +
  scale_y_discrete(breaks = ynames, drop = T)
b
# 
c <- iceberg_w_preds |> 
  dplyr::count(rf) |> 
  ggplot(aes(x = n, y = fct_rev(rf))) +
  geom_col() +
  geom_text(aes(label = n), nudge_x = 800) +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank()) +
  labs(x = 'count') +
  lims(x = c(0, 18000)) +
  scale_y_discrete(breaks = ynames, drop = F)

a+b+c
rm(a,b,c)


# table of prediction classes across models
pred_table <- iceberg_w_preds |> 
  group_by(type, name, id) |> 
  dplyr::count(rf, glmnet, knn) |> 
  ungroup()

# any elements with Rit classes, where models all agree
has_Rit_proteins <- pred_table |> 
  filter(knn == glmnet, knn == rf) |> 
  filter(str_detect(knn, 'Rit') | str_detect(glmnet, 'Rit') |
           str_detect(rf, 'Rit'))

# filter to find elements containing Rits.
results_with_RITs <- 
  iceberg_w_preds |> 
  filter(name %in% has_Rit_proteins$name & id %in% has_Rit_proteins$id) |> 
  filter(knn != 'Other' & glmnet != 'Other') |> 
  rename_with(~paste0(.x, '_pred'), .cols = c(knn, glmnet)) |> 
  select(-contains('seq'), -c(Arch1:Xer), -database,
         -contains('rf'), -contains('pos')) |> 
  group_by(name, id) |> 
  mutate(type = paste(unique(type), collapse = ', ')) |> 
  distinct()

# results_with_RITs |> View()
write_csv(results_with_RITs, './results/iceberg_RITs.csv')

# filter results to find all elements with integrases
iceberg_classes <- 
  iceberg_w_preds |> 
  filter(!(knn == 'Other' & glmnet == 'Other')) |> 
  rename_with(~paste0(.x, '_pred'), .cols = c(knn, glmnet)) |> 
  select(-contains('seq'), -c(Arch1:Xer), -contains('rf')) |> 
  group_by(name, id) |> 
  mutate(type = paste(unique(type), collapse = ', ')) |> 
  distinct()
# 
# iceberg_classes |> View()

write_csv(iceberg_classes, './results/iceberg_all_integrases.csv')



rm(knn_model, glmnet_model, knn_preds, glmnet_preds, knn_prob, glmnet_prob,
   rf_model, rf_preds, rf_prob, ice_scored, iceberg_classes, iceberg_preds,
   iceberg_w_preds, pred_table, fasta, fasta_path, hmm_folder, ynames, temp_dir,
   hmmsearches, ice_proteins)

## PART 2 Examine RIT elements ------

# Nicole wants to know if the Rits in these strains are identical (similarity)...


# reorganize dataset
rits <- 
  results_with_RITs |> 
  left_join(
    iceberg |> 
      select(-type), by = c("name", "id", "accession", "description", 
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
  'ICE-GI3', 'CAP41785.1', 'CAP41786.1', 'CAP41789.1',  'ABC', '',
  'ICE-GI3', 'CAP41836.1', 'CAP41837.1', 'CAP41838.1',  'ABC', 'contains transposases', 
  'ICEMlSym(NZP2037)', 'ANN60838.1', 'ANN60839.1', 'ANN60840.1', 'CBA', '',
  'ICEMlSym(NZP2037)', 'ANN60865.1', 'ANN60866.1', 'ANN60867.1', 'ABC', ''
)

glimpse(elements)

# prep the protein sequences to join to elements' p1-3
proteins <- 
  iceberg |> 
  ungroup() |> 
  filter(prot_accession %in% c(elements$p1_acc, elements$p2_acc, 
                               elements$p3_acc))
p1 <- proteins |> 
  transmute(
    p1_acc = prot_accession,
    p1_seq = prot_seq)
p2 <- proteins |> 
  transmute(
    p2_acc = prot_accession,
    p2_seq = prot_seq,
  )
p3 <- proteins |> 
  transmute(
    p3_acc = prot_accession,
    p3_seq = prot_seq
  )

# join the protein data to elements, join dna data to elements.
rit_elements <- elements |> 
  left_join(p1) |> left_join(p2) |> left_join(p3) |> 
  distinct()  |> 
  left_join(
    rits |> select(contains('ice'), contains('dna')),
    by = 'ice_name') |> 
  distinct()

rm(proteins, p1, p2, p3, elements)

glimpse(rit_elements)

rit_elements |> pull(dna_acc) |> unique()
rit_elements |> pull(dna_length) |> unique()

# translation codes are all 11, add strains
rit_elements <- rit_elements |> 
  left_join(
    tribble(
      ~dna_acc, ~dna_code,
      'AM902716', 11, 
      'CP016079', 11
    )) |> 
  distinct()  |> 
  left_join(rits |> select(strain, ice_id) |> distinct())

# translate dna string for all reading frames, code = 11
make_translations <- function(dna){
  dna <- s2c(dna)
  tibble(
    frame = glue('transl_{1:6}'),
    offset = rep(c(0,1,2), 2), 
    orientation = rep(c('F','R'), each = 3),
    translation = map2_chr(
      .x = offset,
      .y = orientation,
      .f = ~{
        seqinr::translate(dna, frame = .x, sens = .y, numcode = 11) |> 
          c2s()
      }) 
  )
  }

# locate proteins in translations.
locate_protein <- function(translations, protein, tag){
  find_pos <- 
    translations |> 
    mutate(
      prot_acc = tag,
      start = str_locate(translation, protein)[,1],
      end = str_locate(translation, protein)[,2],
      length = end - start + 1,
    )
  return(find_pos)
}


#  map proteins to gene locations
translations <- map(rit_elements$dna_seq, make_translations)
locations <- pmap(
  .l = list(translations,  rit_elements$p1_seq, rit_elements$p1_acc), 
  .f = function(x,y,z) locate_protein(x,y,z)
  )
locations |> bind_rows() |> select(-translation) |> View()
  
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
