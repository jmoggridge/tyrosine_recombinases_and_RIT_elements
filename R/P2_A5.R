## P2_A5 figure out how to identify rit element from genbank...

## Setup ----

library(tidyverse)
library(progress)
library(glue)
library(Biostrings)
library(furrr)
library(tidymodels)

# for parse_genbank function
source('./R/P2_parse_genbank_xml.R')


## Functions for opening genbank and classifying proteins -----

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

extract_features_table <- function(gbk){
  # unnests the features table of genbank record gbk to return CDS items
  gbk |>
    select(feature_table) |> 
    unnest(feature_table) |> 
    select(-c(gb_feature_key, locus_tag, transl_table)) 
}

## rit_finder functions -----

# first, some data that is necessary for classify_proteins()'s hmmsearches
hmmsearch_filler <- read_rds('./data/hmmsearch_filler.rds')


read_hmmsearch <- function(path){
  # read hmmsearch tblout result for one dataset vs one HMM; 
  # return df of (acc, <hmm_name>)
  read_delim(
    file = path, 
    na = '-', 
    delim = ' ', 
    comment = '#',  
    trim_ws = T,
    col_types = cols(), 
    col_names = FALSE
  ) |> 
    transmute(acc = X1, hmm_name = X3, best_dom_score = X9) |> 
    # keep only the best hmm score for each protein
    group_by(acc) |> 
    filter(best_dom_score == max(best_dom_score)) |> 
    ungroup() |>  
    distinct() |> 
    # name the values column after the subfamily hmm
    pivot_wider(names_from = hmm_name, values_from = best_dom_score) 
}

join_hmmsearches2 <- function(df, files){
  # df must have identifiers column to match hmmer output 
  # (need to match the headers in the fasta passed to hmmsearch)
  # read and combine hmmsearch data with reduce
  # join scores to input data & replace any NAs with zeros
  files |> 
    map(read_hmmsearch) |> 
    purrr::reduce(.f = full_join, by = 'acc') |> 
    right_join(df, by = 'acc') |>
    mutate(across(Arch1:Xer, ~replace_na(.x, 0))) |> 
    relocate(everything(), Arch1:Xer)
}


# classify proteins from unnested genbank feature table,
# returns the integrase classification from the model
classify_proteins <- function(ft){
  
  # make temp fasta file of seqs to score against hmm library
  # include filler sequences at head of fasta file
  fasta_path <- tempfile()
  sequences <- c(hmmsearch_filler$seq, ft$translation)
  fasta <- Biostrings::AAStringSet(sequences)
  names(fasta) <- c(hmmsearch_filler$name, ft$protein_id)
  Biostrings::writeXStringSet(fasta, fasta_path)  
  
  # setup paths for hmms and table outputs; 
  # build hmmsearch calls (but not executed until next step)
  temp_dir <- tempdir()
  junk <- tempfile()
  hmmer_call <- 'hmmsearch --noali -o {junk} --tblout'
  hmmsearches <- 
    tibble(hmm_path = Sys.glob(glue('./models/HMM/*'))) |> 
    mutate(
      out_path = hmm_path |> 
        str_replace(glue('./models/HMM/'), glue('{temp_dir}/')) |> 
        str_replace('\\.hmm', '.tbl'),
      calls = glue('{hmmer_call} {out_path} {hmm_path} {fasta_path}')
    )
  
  # do hmmsearch calls in parallel
  plan(multisession, workers = availableCores())
  hmmsearches <- hmmsearches |> 
    mutate(hmmsearches = future_map_dbl(calls, ~ system(.x))) 
  
  # collect hmmsearch results from temp_dir
  ft_scored <- ft |> 
    mutate(acc = protein_id) |> 
    join_hmmsearches2(files = hmmsearches$out_path) |>
    relocate(-c(Arch1:Xer)) |> 
    arrange(feat_id)

  # make integrase prediction with classifiers, find consensus
  ft_classed <- ft_scored |> 
    bind_cols(
      knn_pred = predict(knn_model, new_data = ft_scored) |> 
        pull(.pred_class),
      glmnet_pred = predict(glmnet_model, new_data = ft_scored) |> 
        pull(.pred_class)
    ) |> 
    select(-c(Arch1:Xer)) |> 
    mutate(consensus_pred = case_when(
      knn_pred == glmnet_pred ~ paste(knn_pred),
      TRUE ~ 'no consensus'
    ))
}

edit_feature_table <- function(ft){
  # determine protein orientations and start/stop 
  ft_edit <- ft |> 
    rename_with(~str_remove(.x,'gb_interval_'), contains('gb')) |> 
    mutate(across(c(from, to), ~as.numeric(.x))) |> 
    mutate(
      orientation = ifelse(from < to, 'forward', 'reverse'),
      start = ifelse(orientation == 'forward', from, to),
      stop =  ifelse(orientation == 'forward', to, from)
    ) |> 
    select(-c(from, to, knn_pred, glmnet_pred)) |> 
    mutate(across(c(start, stop), ~as.numeric(.x))) |> 
    relocate(feat_id, everything(), -codon_start)
  return(ft_edit)
}

pivot_feature_table <- function(ft_edit){
  # for each feature (besides the first and last)
  # take info from previous and next features.
  ft_edit |> 
    transmute(
      p1_id = lag(protein_id),
      p1_pred = lag(consensus_pred),
      p1_product = lag(product),
      p1_start = lag(start),
      p1_stop = lag(stop),
      p1_orientation = lag(orientation),
      p2_id = protein_id,
      p2_pred = consensus_pred,
      p2_product = product,
      p2_start = start,
      p2_stop = stop,
      p2_orientation = orientation,
      p3_id = lead(protein_id),
      p3_pred = lead(consensus_pred),
      p3_product = lead(product),
      p3_start = lead(start),
      p3_stop = lead(stop),
      p3_orientation = lead(orientation),
      lflank_id = lag(protein_id, 2),
      lflank_pred = lag(consensus_pred, 2),
      lflank_product = lag(product, 2),
      lflank_start = lag(start, 2),
      lflank_stop = lag(stop, 2),
      lflank_orientation = lag(orientation, 2),
      rflank_id = lead(protein_id, 2),
      rflank_pred = lead(consensus_pred, 2),
      rflank_product = lead(product, 2),
      rflank_start = lead(start, 2),
      rflank_stop = lead(stop, 2),
      rflank_orientation = lead(orientation, 2),
      p1_adjacent = ifelse(feat_id == lag(feat_id) + 1, yes = T, no = F),
      p3_adjacent = ifelse(feat_id == lead(feat_id) - 1, yes = T, no =  F),
    ) |> 
    filter(!is.na(p1_id), !is.na(p3_id)) |> 
    filter(p2_pred != 'Other') 
}


rit_tester <- function(ft_pivot){
  
  # investigate gaps, CDS trio length, whether CDS are Rit subfamilies.
  ft_pivot |> 
    mutate(across(contains('pred'), ~as.character(.x))) |> 
    # check gap sizes, full span of trio CDSs
    mutate(
      # find size of gap/ovelap between CDS
      dist_p1_p2 = p2_start - p1_stop,
      dist_p2_p3 = p3_start - p2_stop,
      # test whether both gaps smaller than 100 bp
      # test length of element
      element_CDS_length = p3_stop - p1_start
    ) |> 
    rowwise() |> 
    # tests
    mutate(
      # both gaps < 250 bp? overlaps are negative.
      rit_dist_check = all(dist_p1_p2 < 250, dist_p2_p3 < 250),
      # full element is less than 4kbp?
      rit_length_check = ifelse(element_CDS_length < 4000, T, F),
      # all CDS are predicted Rit subfamilies members
      rit_all_check = all(
        map_lgl(c(p1_pred, p2_pred, p3_pred), ~str_detect(.x, 'Rit'))
      ),
      # trio has ABC or CBA arrangement?
      rit_ABC_check = any(
        all(p1_pred == 'RitA', p2_pred == 'RitB', p3_pred == 'RitC'),
        all(p1_pred == 'RitC', p2_pred == 'RitB', p3_pred == 'RitA')
      ),
      # trio shares orientation?
      rit_orientation = case_when(
        all(p1_orientation == 'forward', 
            p2_orientation == 'forward',
            p3_orientation == 'forward') ~ 'all forward',
        all(p1_orientation == 'reverse', 
            p2_orientation == 'reverse',
            p3_orientation == 'reverse') ~ 'all reverse',
        TRUE ~ 'not same orientation'
      )
    )
}

rit_selector <- function(rit_tests, ft_edit){
  # rit selector: check tests, filter candidates,
  # package up upstream and downstream cds.
  rits <- rit_tests |> 
    rowwise() |> 
    filter(all(rit_dist_check, rit_all_check, rit_length_check)) |> 
    mutate(
      upstream_cds = map(
        p1_start, 
        ~ft_edit |> 
          filter(.x - stop > -200 & .x - stop < 1000)
      ),
      downstream_cds = map(
        p3_stop, 
        ~ft_edit |> 
          filter(start -.x > -200 & start - .x < 1000)
      )
    ) |> 
    ungroup() |> 
    nest(rits = everything()) |> 
    transmute(nuc_id = x, rits)
  return(rits)
}

## Main ----

## load classifiers
knn_model <- read_rds('./models/knn_classifier.rds')
glmnet_model <- read_rds('./models/glmnet_classifier.rds')

### Load probable RIT elements id data ----

# genbank files for these nuc ids lack the CDS features.
no_cds <- c('1834417989', '1024364864')

# these nucleotides probably contain RIT element bc have 3 CDD proteins.
three_ints <- read_rds('./data/CDD/ids_w_three_integrases.rds') |> 
  filter(!nuc_id %in% no_cds) |> 
  arrange(slen)
glimpse(three_ints)

# prot_id and accessions for this set of nt's.
three_ints_prots <- three_ints |> unnest(prot_id) |> pull(prot_id)

# join with other accession numbers
prot_accessions <- read_rds('./data/CDD/prot_summary_fixed.rds') |> 
  filter(prot_id %in% three_ints_prots) |> 
  transmute(prot_id, prot_accession = glue('{caption}.1'))

# join the protein accession labels that match the genbank CDS entries
three_ints <- three_ints |> 
  unnest(c(prot_id, cdd_title)) |> 
  left_join(prot_accessions, by = c("prot_id")) |> 
  group_by(nuc_id) |> 
  nest(prot_data = c(prot_id, prot_accession, cdd_title)) |> 
  mutate(nuc_accession = caption)

# create an index for nuc ids to their genbank records
genbank_files_index <- 
  tibble(file = Sys.glob('./data/CDD/RIT_gbk_[0-9].rds')) |> 
  mutate(nuc_id = map(file, ~ read_rds(.x) |> pull(nuc_id))) |> 
  unnest(nuc_id)

rm(prot_accessions, three_ints_prots, no_cds)

##  Find RITs -----


# find RIT elements for nucleotide id x.
# x's genbank file needs to be in genbank_files_index for # open_genbank.
# CDS extracted by extract_features. 
# Protein CDSs are classified by classify_proteins. 
# Then classified proteins scanned for RIT arrangements

rit_finder <- function(x){
  
  # open genbank and extract features table for nuc_id x,
  ft <- open_genbank(x = x) |> 
    extract_features_table() 
  
  # predict integrase classes for all protein CDSs
  ft_class <- ft |>
    filter(!is.na(translation)) |> 
    classify_proteins() |>
    select(-acc)
  # join predictions to feature table
  ft_join <- ft |> 
    left_join(ft_class,
              by = c("gb_interval_from", "gb_interval_to", "feat_id", 
                     "codon_start", "product", "protein_id", "translation")) |> 
    mutate(pseudo_or_partial = ifelse(is.na(protein_id), T, F))
  
  # determine protein orientations and start/stop; drop some data 
  ft_edit <- edit_feature_table(ft = ft_join)
  
  # for each feature (besides the first and last)
  # take info from previous and next features.
  ft_pivot <- pivot_feature_table(ft_edit = ft_edit)
  rit_tests <- rit_tester(ft_pivot = ft_pivot)
  
  # rit selector: check tests, filter candidates,
  # package up upstream and downstream cds.
  rits <- rit_selector(rit_tests = rit_tests, ft_edit = ft_edit)
  return(rits)
  }

# example

x <- three_ints$nuc_id[[2]]
x
df <- rit_finder(x)
df


# TODO figure out the issue with nuc_id 2
## issue was HMM search results were empty so had no HMM name in col X3


# are the three RIT proteins there??
# prot_x <- three_ints |>
#   filter(nuc_id == x) |> 
#   unnest(prot_data) |> 
#   pull(prot_accession)
# 
# ft[which(ft$protein_id %in% prot_x), ] |> 
#   left_join()



# figure out where there is three RITs in a row...



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
