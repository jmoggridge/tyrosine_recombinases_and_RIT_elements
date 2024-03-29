## P2  rit_finder algorithm

library(tidyverse)
library(progress)
library(glue)
library(Biostrings)
library(furrr)
library(tidymodels)

# load parse_genbank function
source('./R/P2_parse_genbank_xml.R')

# classifiers
knn_model <- read_rds('./models/knn_classifier.rds')
glmnet_model <- read_rds('./models/glmnet_classifier.rds')

# some data that is necessary for classify_proteins()'s hmmsearches
hmmsearch_filler <- read_rds('./data/hmmsearch_filler.rds')


## rit_finder functions -----


# for a nuc_id, return the parsed genbank record as a tibble
# genbank_files_index needs to be loaded in the environment

# open_genbank <- function(x){
#   genbank_files_index |> 
#     # figure out which file to open from index
#     filter(nuc_id == x) |> 
#     pull(file) |> 
#     # pull the genbank record for id and parse it
#     read_rds() |> 
#     filter(nuc_id == x) |> 
#     unnest(gbk) |> 
#     pull(gbk) |> 
#     parse_genbank()
# }

open_genbank <- function(x, genbank_lib){
  # pull the genbank record for id and parse it
  read_rds(genbank_lib) |> 
    filter(nuc_id == x) |> 
    unnest(gbk) |> 
    pull(gbk) |> 
    parse_genbank()
}

## TODO secondary parser for downloaded genbank .xml files
# open_genbank2 <- function(x){
#   gb <- read_file(glue('./data/CDD/genbank_w_cds/{x}.xml'))
#   gb |> parse_genbank()
# }


extract_features_table <- function(ft){
  # unnests the features table of genbank record gbk to return CDS items
  ft |> select(feature_table) |> 
    unnest(feature_table) |> 
    select(-c(matches('gb_feature_key'), 
              matches('locus_tag'),
              matches('transl_table')) 
           )
}

read_hmmsearch <- function(path){
  # read hmmsearch tblout result for one dataset vs one HMM; 
  # return df of (acc, <hmm_name>)
  read_delim(
    file = path, na = '-', delim = ' ', comment = '#',  trim_ws = T,
    col_types = cols(), col_names = FALSE
  ) |> 
    transmute(acc = X1, hmm_name = X3, best_dom_score = X9) |> 
    # keep only the best hmm score for each protein
    group_by(acc) |> filter(best_dom_score == max(best_dom_score)) |> 
    ungroup() |>  distinct() |> 
    # name the values column after the subfamily hmm
    pivot_wider(names_from = hmm_name, values_from = best_dom_score) 
}

join_hmmsearches2 <- function(df, files){
  # df must have identifiers column to match hmmer output 
  # (need to match the headers in the fasta passed to hmmsearch)
  # read and combine hmmsearch data with reduce
  # join scores to input data & replace any NAs with zeros
  files |> map(read_hmmsearch) |> 
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
  hmmer_call <- glue('hmmsearch --noali -o {junk} --tblout')
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
  if(!'product' %in% colnames(ft_edit)){
    ft_edit <- ft_edit |> 
      mutate(product = as.character(NA))
  }
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
      # proteins flanking upstream and downstream
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
      # are the p1, p2, and p3 CDSs in series (T) or seperated by other genes(F)?
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
      # TODO FIX OVERLAP
      # find size of gap/ovelap between CDS
      overlap_p1p2 =  p1_stop - p2_start + 1,
      overlap_p2p3 = p2_stop - p3_start + 1,
      # test whether both gaps smaller than 100 bp
      # test length of element
      element_CDS_length = p3_stop - p1_start
    ) |> 
    rowwise() |> 
    # tests
    mutate(
      # both gaps < 250 bp? overlaps are negative.
      rit_dist_check = all(overlap_p1p2 > -250, overlap_p2p3 > -250), 
      rit_overlap_check = all(overlap_p1p2 < 200, overlap_p2p3 < 200), 
      # full element is less than 4kbp?
      rit_length_check = ifelse(element_CDS_length < 8000, T, F),
      # CDS trio are all predicted integrases 
      rit_integrase_check = all(
        map_lgl(c(p1_pred, p2_pred, p3_pred), ~!str_detect(.x, 'Other'))
      ),      
      # CDS trio are all predicted Rit subfamilies members
      rit_all_check = all(
        map_lgl(c(p1_pred, p2_pred, p3_pred), ~str_detect(.x, 'Rit'))
      ),
      # CDS trio has RIT_ABC or CBA arrangement?
      rit_ABC_check = any(
        all(p1_pred == 'RitA', p2_pred == 'RitB', p3_pred == 'RitC'),
        all(p1_pred == 'RitC', p2_pred == 'RitB', p3_pred == 'RitA')
      ),
      # CDS trio shares gene orientation?
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


# rit selector: check tests, filter candidates,
rit_selector <- function(nuc_id, rit_tests, ft_edit){
  rits <- rit_tests |> 
    rowwise() |> 
    # keep trios that pass: all - integrases, 
    filter(all(rit_dist_check, rit_overlap_check,
               rit_integrase_check, rit_length_check)) |> 
    # package up upstream and downstream cds.
    mutate(
      upstream_cds = map(
        p1_start, 
        ~ft_edit |> filter(.x - stop > -200 & .x - stop < 1000)
      ),
      downstream_cds = map(
        p3_stop, 
        ~ft_edit |> filter(start -.x > -200 & start - .x < 1000)
      )
    ) |> 
    ungroup() |> 
    nest(rits = everything()) |> 
    transmute(rits, success = T)
  return(rits)
}


##  Find RITs -----


# find RIT elements for nucleotide id x.
# x's genbank file needs to be in genbank_files_index for # open_genbank.
# CDS extracted by extract_features. 
# Protein CDSs are classified by classify_proteins. 
# Then classified proteins scanned for RIT arrangements

rit_finder <- function(x, genbank_lib){
  
  # open genbank and extract features table for nuc_id x,
  gbk <- open_genbank(x = x, genbank_lib = genbank_lib)
  
  # if no CDS in feature table, return NA
  # TODO use 2ndary parser for downloaded .xml files
  
  if(any(class(gbk) == 'logical')) {
    ## try alternate parser??
    # ft <- open_genbank2(x)
    return(tibble(nuc_id_fail = x, rits = list(NA), success = F))
  }
  # extract CDS features
  ft <- extract_features_table(gbk) 
  
  # predict integrase classes for all protein CDSs
  ft_class <- ft |>
    filter(!is.na(translation)) |> 
    classify_proteins() |>
    select(-acc)
  
  # join predictions to feature table
  ft_join <- ft |> 
    left_join(ft_class) |> 
    mutate(pseudo_or_partial = ifelse(is.na(protein_id), T, F))
  
  # determine protein orientations and start/stop; drop some data 
  ft_edit <- edit_feature_table(ft = ft_join)
  
  # for each feature (besides the first and last)
  # take info from previous and next features.
  ft_pivot <- pivot_feature_table(ft_edit = ft_edit)
  rit_tests <- rit_tester(ft_pivot = ft_pivot)
  
  # rit selector: check tests, filter candidates,
  # package up upstream and downstream cds.
  rits <- rit_selector(
    rit_tests = rit_tests, 
    ft_edit = ft_edit, 
    nuc_id = x
    )
  
  return(rits)
}


# same as rit funder one but directly from the genbank containing tibble
# parse genbank xml string gbk_xml and extract features table for nuc_id x
rit_finder2 <- function(x, gbk_xml){
  
  # parse genbank xml string to tibble
  gbk <- parse_genbank(gbk_xml)
  
  # if no CDS in feature table, return NA
  # TODO use 2ndary parser for downloaded .xml files
  if(any(class(gbk) == 'logical')) {
    ## try alternate parser??
    # ft <- open_genbank2(x)
    return(tibble(nuc_id_fail = x, rits = list(NA), success = F))
  }
  # extract CDS features
  ft <- extract_features_table(gbk) 
  
  # predict integrase classes for all protein CDSs
  ft_class <- ft |>
    filter(!is.na(translation)) |> 
    classify_proteins() |>
    select(-acc)
  
  # join predictions to feature table
  ft_join <- ft |> 
    left_join(ft_class) |> 
    mutate(pseudo_or_partial = ifelse(is.na(protein_id), T, F))
  
  # determine protein orientations and start/stop; drop some data 
  ft_edit <- edit_feature_table(ft = ft_join)
  
  # for each feature (besides the first and last)
  # take info from previous and next features.
  ft_pivot <- pivot_feature_table(ft_edit = ft_edit)
  rit_tests <- rit_tester(ft_pivot = ft_pivot)
  
  # rit selector: check tests, filter candidates,
  # package up upstream and downstream cds.
  rits <- rit_selector(
    rit_tests = rit_tests, 
    ft_edit = ft_edit, 
    nuc_id = x
  )
  
  return(rits)
}


# rm(gbk, ft, ft_class, ft_edit, ft_pivot, ft_join, rits, rit_tests)
