## Classifier training / assessment by repeated nested-stratified CV

## From training data do k-fold cv
library(tidyverse)
library(tidymodels)
library(Biostrings)
library(furrr)
library(DECIPHER)
library(glue)
library(here)
library(tictoc)

# create directory structure for classifier files (alignments, hmms, results for each resample)
dir <- 'classifier_01/'
out_path <- glue(here::here(), '/data/', dir)
system(glue('mkdir {out_path}'))

list('align','hmm', 'hmmsearch', 'results') |> 
  map(~glue(dir, .x)) |> 
  map(~glue(here::here(), '/data/', .x)) |>
  map(~system(glue('mkdir {.x}')))

rm(dir)

## Combine datasets  --------------------------------------------

# integrase dataset
smart_df <- 
  read_rds('./data/SMART/smart_df.rds') |>
  select(subfamily, acc, description, prot_seq, dom_seq)

# 20 subfamilies for classification
subfamilies <- unique(smart_df$subfamily)

# negative examples (non-integrases)
non_integrases <- 
  read_rds('./data/non_integrase_seqs/nonint_df.rds') |> 
  mutate(dom_seq = NA) |> 
  select(subfamily, acc, description, prot_seq, dom_seq)

full_dataset <- bind_rows(smart_df, non_integrases)
rm(smart_df, non_integrases)

## Initial split -------------------------------------------------
set.seed(54321)

# TODO (remove downsampling)
df <- full_dataset |> 
  group_by(subfamily) |> 
  slice_sample(n = 2000, replace = F) |> 
  ungroup()
df |> count(subfamily) |> print.AsIs()

df_split <- initial_split(df, 0.75, strata = subfamily)

train <- training(df_split)
test <- testing(df_split)

train |> count(subfamily) |> print.AsIs()
test |> count(subfamily) |> print.AsIs()

rm(df, full_dataset)

## full data split
# df_split <- initial_split(full_dataset, 0.75, strata = subfamily)
# train <- training(df_split)
# test <- testing(df_split)
# write_rds(test, glue('{out_path}/results/final_test.rds'))


#### Nested CV --------------------------------------------------

set.seed(1234)
nest_cv <- 
  nested_cv(
    train, 
    outside = vfold_cv(v = 3, repeats = 1, strata = subfamily),
    inside = vfold_cv(v = 3, repeats = 1, strata = subfamily)
  ) |> 
  transmute(
    outer_id = paste0('outer', str_remove(id, 'Fold')),
    outer_splits = splits,
    inner_resamples
  ) |>
  unnest(inner_resamples) |> 
  transmute(
    outer_id, outer_splits,
    inner_id = paste0(outer_id, '_', id),
    inner_splits = splits
  ) |> 
  nest(inner_resamples = c(inner_id, inner_splits))

nest_cv


#### Training Functions ---------------------------------------


## Prepare domains for alignment/hmm 
prep_domains_df <- function(training){
  training |> 
    # subset domain subfamilies
    filter(subfamily !=  'non_integrase') |> 
    select(subfamily, acc, dom_seq) |> 
    # remove any duplicated domain sequences
    group_by(dom_seq) |> sample_n(1) |> ungroup() |> 
    # TODO (remove) downsample subfamilies
    group_by(subfamily) |> 
    slice_max(n = 100, order_by = row_number()) |> 
    ungroup() |> 
    # duplicate subfamily so that info gets passed to map later 
    mutate(subfam = subfamily) |> 
    # nest dataframes by subfamily 
    group_by(subfamily) |> nest() |> ungroup() |> 
    # arrange largest to smallest
    mutate(nrow = map_int(data, nrow)) |> 
    arrange(desc(nrow)) |> 
    select(-nrow) 
}

## Alignment
# pass a dataframe with subfamily, acc, dom_seq
# return same df with alignment stringset column
align_domains <- function(df, dest){
  outpath <- paste0(dest, df$subfam[1], '.aln')
  aa_set <- Biostrings::AAStringSet(df |> pull(dom_seq))
  names(aa_set) <- df$acc
  aligned <- DECIPHER::AlignSeqs(aa_set, verbose = F, processors = 1)
  Biostrings::writeXStringSet(aligned, filepath = outpath)
  return(aligned)
}

# Wrapper to apply align_domains to each subfamily
build_alignments_library <- function(domains, fold){
  path <- glue('{out_path}/align/', fold, '/')
  # do align_domains in parallel
  plan(multisession, workers = availableCores() - 1)
  domains |> mutate(
      aligned = future_map(
        .x = data, 
        .f = ~ align_domains(.x, dest = path),
        .options = furrr_options(scheduling = Inf)
      )
    )
}
  

# create paths, build call, do calls to hmmbuild for each subfamily
build_hmm_library <- function(fold){
  tibble(
      align_in_path = Sys.glob(glue('{out_path}align/', fold, '/*'))
    ) |> 
    mutate(
      hmm_out_path = str_replace_all(align_in_path, 'align|aln', 'hmm'),
      hmmbuild_call = glue('hmmbuild {hmm_out_path} {align_in_path}'),
      hmmbuild = map_dbl(hmmbuild_call, ~ system(.x))
    )
  }

# score a dataframe of sequences against HMM library with hmmsearch
hmmsearch_scores <- function(df, fold){
  
  # make temp fasta file of seqs to score against hmm library
  temp <- tempfile()
  fasta <- Biostrings::AAStringSet(df$prot_seq)
  names(fasta) <- df$acc
  writeXStringSet(fasta, temp)  
  
  # setup paths for hmms and table outputs; build hmmsearch calls
  temp <- tempfile()
  hmmsearches <- 
    tibble(hmm_path = Sys.glob(glue('{out_path}hmm/', fold, '/*'))) |> 
    mutate(
      out_path = hmm_path |> 
        str_replace('/hmm/', '/hmmsearch/') |> 
        str_replace('\\.hmm', '.tbl'),
      files = glue('{out_path} {hmm_path} {temp}'),
      calls = glue('hmmsearch --noali -o {temp} --tblout {files}')
    )
  plan(multisession, workers = availableCores())
  hmmsearches <- hmmsearches |> 
    mutate(hmmsearches = future_map_dbl(calls, ~ system(.x))) 
  return(hmmsearches)
}



inner_fit <- function(fold){
  
  system(glue('mkdir {out_path}/align/', fold))
  system(glue('mkdir {out_path}/hmm/', fold))
  system(glue('mkdir {out_path}/hmmsearch/', fold))
  
  # prep
  training <- analysis(inner_fold)
  testing <- assessment(inner_fold)
  domains <- 
    analysis(inner_fold)
    prep_domains_df()
  
  # align
  message('Aligning training domains...')
  tic()
  aligns <- build_alignments_library(domains, fold)
  toc()
  
  # train hmm
  message('Building HMMs...')
  hmm_check <- build_hmm_library(fold)
  message('hmmbuild calls all finished without error?')
  all(hmm_check$hmmbuild == 0)
  
  # hmm scores test and train
  message('Scoring sequences against HMM library...')
  train_search <- hmmsearch_scores(training, fold)
  message('hmmsearch calls all finished without error?')
  all(train_search$hmmsearches == 0)
  
  
  # train classifier
  # evaluate classifier
}


### MAIN / TRAIN & TEST -----
rm(train, test, df_split)
inner_cv <- nest_cv$inner_resamples[[2]]
inner_cv
inner_fold <- inner_cv$inner_splits[[2]]
inner_fold
fold <- inner_cv$inner_id[[2]]
fold

# make directories for this fold's alignment and hmm libraries
system(glue('mkdir {out_path}/align/', fold))
system(glue('mkdir {out_path}/hmm/', fold))
system(glue('mkdir {out_path}/hmmsearch/', fold))

training <- analysis(inner_fold)
testing <- assessment(inner_fold)
domains <- inner_fold |> analysis() |> prep_domains_df()

# check stratification is adequate / no empty classes
message("Train and test split contain > 10 obs per subfamily?")
all(
  analysis(inner_fold) |> count(subfamily) |> pull(n) %>% all(.data > 10), 
  assessment(inner_fold) |> count(subfamily) |> pull(n) %>% all(.data > 10)
)

tic()
aligns <- build_alignments_library(domains, fold)
toc()

hmm_check <- build_hmm_library(fold)
all(hmm_check$hmmbuild == 0)

train_search <- hmmsearch_scores(training, fold)
all(train_search$hmmsearches == 0)

# parse hmmscores

## TODO finish hmmscoring & parsing

# parse HMM search tables
hmmsearches |> 
  mutate(hmmsearch_tbl = )




##-----





# make_hmms

# glob <- 
#   tibble(
#     align_in_path = Sys.glob(glue('{out_path}align/', fold, '/*'))
#   ) |> 
#   mutate(
#     hmm_out_path = str_replace_all(align_in_path, 'align|aln', 'hmm'),
#     hmmbuild_call = glue('hmmbuild {hmm_out_path} {align_in_path}')
#   )
#   
# glob |> mutate(hmmbuild = map_dbl(hmmbuild_call, ~ system(.x)))





### HMMSEARCH SCORING 
# df <- training
# 
# # make temp fasta file of seqs to score against hmm library
# temp <- tempfile()
# fasta <- Biostrings::AAStringSet(df$prot_seq)
# names(fasta) <- df$acc
# writeXStringSet(fasta, temp)  
# 
# # setup paths for hmms and table outputs; build hmmsearch calls
# hmmsearches <- 
#   tibble(hmm_path = Sys.glob(glue('{out_path}hmm/', fold, '/*'))) |> 
#   mutate(
#     out_path = hmm_path |> 
#       str_replace('/hmm/', '/hmmsearch/') |> 
#       str_replace('\\.hmm', '.tbl'),
#     files = glue('{out_path} {hmm_path} {temp}'),
#     calls = glue('hmmsearch --noali -o temp --tblout {files}')
#   )
# 
# plan(multisession, workers = availableCores())
# hmmsearches <- hmmsearches |> 
#   mutate(hmmsearches = future_map_dbl(hmmsearch_call, ~ system(.x))) 
