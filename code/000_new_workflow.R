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

##  Plots ------------------------------------------------------------

# plot scores 
plt_hmmscores_boxplots <- function(df, title){
  df |> 
    pivot_longer(Arch1:Xer) |> 
    mutate(fam = ifelse(subfamily == name, T, F)) |> 
    ggplot(aes(value, fct_rev(name), color = fam)) +
    geom_boxplot(outlier.alpha = 0.3, outlier.size = 0.5, outlier.shape = 0) +
    facet_wrap(~subfamily, nrow = 4) +
    labs(y ='', x = 'HMM score', color = '',
         subtitle = glue('{title}'),
         caption = 'Each subfamily is shown as a panel with the distributions \
         of scores for each HMM.') +
    theme_light()
}
# plt_hmmscores_boxplots(training, 'Training fold')  
# plt_hmmscores_boxplots(testing, 'Testing fold')  


#### Training Functions ---------------------------------------

## Prepare domains for alignment/hmm 
prep_domains_df <- function(training){
  training |> 
    # subset domain subfamilies
    filter(subfamily !=  'non_integrase') |> 
    select(subfamily, acc, dom_seq) |> 
    # remove any duplicated domain sequences
    group_by(dom_seq) |> sample_n(1) |> ungroup() |> 
    # TODO ****(remove) downsample subfamilies****
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

# score a dataframe of sequences against the HMM library with hmmsearch
hmmsearch_scores <- function(df, fold, tag){
  
  # make temp fasta file of seqs to score against hmm library
  temp <- tempfile()
  fasta <- Biostrings::AAStringSet(df$prot_seq)
  names(fasta) <- df$acc
  writeXStringSet(fasta, temp)  
  
  # setup paths for hmms and table outputs; build hmmsearch calls
  hmmsearches <- 
    tibble(hmm_path = Sys.glob(glue('{out_path}hmm/', fold, '/*'))) |> 
    mutate(
      out_path = hmm_path |> 
        str_replace('/hmm/', '/hmmsearch/') |> 
        str_replace('\\.hmm', glue('.{tag}.tbl')),
      calls = glue(
        'hmmsearch --noali -o temp --tblout {out_path} {hmm_path} {temp}'
        )
    )
  # do hmmsearch calls in parallel
  plan(multisession, workers = availableCores())
  hmmsearches <- hmmsearches |> 
    mutate(hmmsearches = future_map_dbl(calls, ~ system(.x))) 
  return(hmmsearches)
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
    group_by(acc) |> 
    filter(best_dom_score == max(best_dom_score)) |> 
    ungroup() |>  
    distinct() |> 
    # name the values column after the subfamily hmm
    pivot_wider(names_from = hmm_name, values_from = best_dom_score) |> 
    unnest(cols = c())
}

# reads all hmmsearch tables in glob and joins data
join_hmmsearches <- function(df, files){
  searches <- 
    # read and combine hmmsearch data
    map(files, read_hmmsearch) |> 
    purrr::reduce(.f = full_join, by = 'acc') 
  # join scores to input data & replace any NAs with zeros
  df <- df |> 
    left_join(searches, by = 'acc')
  df <- df |>
    mutate(across(Arch1:Xer, ~replace_na(.x, 0))) 
}




## Wrappers for CV ----------------


#' inner_prep: function to do an inner cv fold training
#' calls prep_domains_df, build_alignments_library, build_hmm_library,
#' hmmsearch_scores, and join_hmmsearches
inner_prep <- function(split, id){
  
  print(glue("\n\nPreparing fold: {id} --------- \n"))
  print(split)
  
  system(glue('mkdir {out_path}/align/', id))
  system(glue('mkdir {out_path}/hmm/', id))
  system(glue('mkdir {out_path}/hmmsearch/', id))
  
  # check stratification is adequate / no empty classes
  message("\n\nTrain and test split contain > 10 obs per subfamily?")
  all(
    analysis(split) |> count(subfamily) |> pull(n) %>% all(.data > 10), 
      assessment(split) |> count(subfamily) |> pull(n) %>% all(.data > 10)
    ) |> 
    print()
  
  # align
  message('\n\nAligning training domains...')
  tic()
  aligns <- 
    analysis(split) |> 
    prep_domains_df() |> 
    build_alignments_library(fold = id)
  toc()
  
  # train hmm
  message('\n\nBuilding HMMs...')
  hmm_check <- build_hmm_library(id)
  message('\n\nhmmbuild calls all finished without error?')
  print(all(hmm_check$hmmbuild == 0))
  
  # hmm scores test and train
  message('\n\nScoring sequences against HMM library...')
  tic()
  train_search <- 
    analysis(split) |> 
    hmmsearch_scores(fold = id, tag = 'train')
  test_search <- 
    assessment(split) |> 
    hmmsearch_scores(fold = id, tag = 'test')
  toc()
  
  message('\n\nHMM searches all finished without error?')
  all(train_search$hmmsearches == 0)
  
  # parse hmmscores and add to data frame
  train <-  
    analysis(split) |> 
    join_hmmsearches(files = train_search$out_path)
  test <-  
    assessment(split) |> 
    join_hmmsearches(files = test_search$out_path)
  
  # return tibble row with train and test data for classifiers
  message('\n\nReturning datasets --------------------------------------')
  tibble(train = list(train), test = list(test))
}

# inner_cv <- function(){}
# outer_cv <- function(outer_split, inner_cv, fold){}
# outer_fit <- 

### MAIN / TRAIN & TEST ----------------------------------------

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
  select(subfamily, acc, description, prot_seq, dom_seq) |> 
  mutate(subfamily = as_factor(subfamily))

# 20 subfamilies for classification
subfamilies <- unique(smart_df$subfamily)

# negative examples (non-integrases)
non_integrases <- 
  read_rds('./data/non_integrase_seqs/nonint_df.rds') |> 
  mutate(dom_seq = NA,
         subfamily = as_factor(subfamily)) |> 
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
    outer_id, 
    outer_splits,
    inner_id = paste0(outer_id, '_', id),
    inner_splits = splits
  ) |> 
  nest(inner_resamples = c(inner_id, inner_splits))

nest_cv |> unnest(inner_resamples)

rm(train, test, df_split)


## DEV inner cv fold --------------------------------------------------

cv <- nest_cv$inner_resamples[[2]]
cv
split <- inner_cv$inner_splits[[2]]
split
id <- inner_cv$inner_id[[2]]
id



# inner_prep(split = split,  id = id)

## TODO modelling....
## TODO setup cv functions


inner_cv <- function(resamples){
  # prepare scored datasets
  resamples |> 
    mutate(inner_prep = map2(inner_splits, inner_id, ~inner_prep(.x, .y))) |> 
    unnest(inner_prep)
  # fit model
  # predict test
  # evaluate performance
  # return model, params, performance mean & 
}


df <- cv |> 
  mutate(inner_prep = map2(inner_splits, inner_id, ~inner_prep(.x, .y)))
df <- df |> unnest(inner_prep)


train <-  df$train[[1]] |> mutate(subfamily = as_factor(subfamily)) 
test <- df$test[[1]] |> mutate(subfamily = as_factor(subfamily)) 


evaluate_model <- function(model, train, test){
  recip <- 
    recipe(subfamily ~ ., data = train) |> 
    update_role(c('acc','description', contains('seq')), new_role = "id") |> 
    step_scale(all_predictors())
  fit_wf <- workflow() |>
    add_model(model) |> 
    add_recipe(recip) |>
    fit(data = train)
  predictions <- fit_wf |> 
    predict(test) %>% 
    bind_cols(test)
  bind_rows(
    mcc(predictions, truth = subfamily, estimate = .pred_class),
    spec(predictions, truth = subfamily, estimate = .pred_class),
    sens(predictions, truth = subfamily, estimate = .pred_class)
  )
}
  
## unfit models --------------------------------------------

# update model specifications using parameters grid
make_models <- function(model, grid){
  models <- map(
    seq_len(nrow(grid)), ~{
      mods <- model |> 
        update(grid[.x, ])
    })
  tibble(spec = models) |> 
    mutate(model_id = row_number()) |> 
    bind_cols(grid)
}

# decision tree models
tree_models <- 
  decision_tree(mode = 'classification') |> 
  set_engine('rpart') |> 
  make_models(grid = grid_max_entropy(
    size = 20, 
    tree_depth(), 
    min_n(range = c(2L, 40L)), 
    cost_complexity()
  ))

# logistic reg models
log_reg_models <- 
  logistic_reg(mode = 'classification') |> 
  set_engine('glmnet') |> 
  make_models(grid = grid_max_entropy(
    mixture(), 
    penalty(), 
    size = 20
  ))
## unfit models --------------------------------------------


nest_cv$inner_resamples[[1]] |> inner_cv()


dt |> mutate(eval = map(model, ~evaluate_model(.x, train, test)))


# 
# # make directories for this fold's alignment and hmm libraries
# system(glue('mkdir {out_path}/align/', fold))
# system(glue('mkdir {out_path}/hmm/', fold))
# system(glue('mkdir {out_path}/hmmsearch/', fold))
# 
# training <- analysis(inner_fold)
# testing <- assessment(inner_fold)
# domains <- inner_fold |> analysis() |> prep_domains_df()
# 
# # check stratification is adequate / no empty classes
# message("Train and test split contain > 10 obs per subfamily?")
# all(
#   analysis(inner_fold) |> count(subfamily) |> pull(n) %>% all(.data > 10), 
#   assessment(inner_fold) |> count(subfamily) |> pull(n) %>% all(.data > 10)
# )
# 
# tic()
# aligns <- build_alignments_library(domains, fold)
# toc()
# 
# hmm_check <- build_hmm_library(fold)
# all(hmm_check$hmmbuild == 0)
# 
# 
# train_search <- hmmsearch_scores(training, fold, tag = 'train')
# all(train_search$hmmsearches == 0)
# test_search <- hmmsearch_scores(testing, fold, tag = 'test')
# all(test_search$hmmsearches == 0)
# 
# 
# # set gathering threshold for each model as the lowest score of true positives
# fit_thresholds <- function(training){
#   train |> 
#     select(subfamily, Arch1:Xer) |> 
#     pivot_longer(cols = Arch1:Xer, 
#                  names_to = 'hmm_name', 
#                  values_to = 'threshold') |> 
#     filter(subfamily == hmm_name) |> 
#     group_by(subfamily) |> 
#     filter(threshold == min(threshold)) |> 
#     ungroup() |> 
#     distinct() |> 
#     select(-subfamily)
# }
# 
# 
# # parse hmmscores and add to data frame
# training <- analysis(inner_fold) |> 
#   join_hmmsearches(files = train_search$out_path)
# testing <- assessment(inner_fold) |> 
#   join_hmmsearches(files = test_search$out_path)
# 
#   # parse hmmscores and add to data frame
# tibble(
#   training = list(analysis(inner_fold) |> 
#     join_hmmsearches(files = train_search$out_path)),
#   testing = list(assessment(inner_fold) |> 
#     join_hmmsearches(files = test_search$out_path))
# )
# 
# training
# testing


# .options = future_options(packages = "parsnip")