## Classifier training / assessment by repeated nested-stratified CV

# Performs a nested CV where the training data are used to create domain subfamily alignments and train HMMs. Then the HMMs aer used to score train and test sequences. The scores are used to predict the class of the test data, and performance metrics for the different models are compared in the inner CV (model tuning) and outer CV (model selection).

# TODO selection of models
# TODO selection of k and rep for outer and inner CV
# TODO finish inner cv function
# TODO write outer cv function

## Libraries ---------------------------------------------------------------
library(tidyverse)
library(tidymodels)
library(Biostrings)
library(furrr)
library(DECIPHER)
library(glue)
library(here)
library(tictoc)
library(crayon)

## Directories ---------------------------------------------------------------

# create directory structure for classifier files (alignments, hmms, results for each resample)
dir <- 'classifier_01/'
out_path <- glue(here::here(), '/', dir)
system(glue('mkdir {out_path}'))

list('align','hmm', 'hmmsearch', 'results') |> 
  map(~glue(dir, .x)) |> 
  map(~glue(here::here(), '/', .x)) |>
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
# TODO (remove downsampling)

set.seed(54321)
df <- full_dataset |> 
  group_by(subfamily) |> 
  slice_sample(n = 500, replace = F) |> 
  ungroup()

df |> count(subfamily) |> print(n=22)

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


## Nested CV --------------------------------------------------

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


# Functions ----

###  Plot functions ------------------------------------------------------------

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


### Prep Functions ---------------------------------------

## Prepare domains for alignment/hmm 
prep_domains_df <- function(training){
  training |> 
    # subset domain subfamilies
    filter(subfamily !=  'non_integrase') |> 
    select(subfamily, acc, dom_seq) |> 
    # remove any duplicated domain sequences
    group_by(dom_seq) |> sample_n(1) |> ungroup() |> 
    
    # TODO ****(remove) downsample subfamilies****
    # group_by(subfamily) |> 
    # slice_max(n = 500, order_by = row_number()) |> 
    # ungroup() |> 
    
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
  
  temp <- tempfile()
  plan(multisession, workers = availableCores())
  
  tibble(msa_path = Sys.glob(glue('{out_path}align/', fold, '/*'))) |> 
    mutate(
      hmm_path = str_replace_all(msa_path, 'align|aln', 'hmm'),
      hmmbuild_call = glue('hmmbuild -o {temp} {hmm_path} {msa_path}'),
      hmmbuild = future_map_dbl(hmmbuild_call, ~ system(.x))
    )
  }

# score a dataframe of sequences against the HMM library with hmmsearch
hmmsearch_scores <- function(df, fold, tag){
  
  # make temp fasta file of seqs to score against hmm library
  fasta_path <- tempfile()
  fasta <- Biostrings::AAStringSet(df$prot_seq)
  names(fasta) <- df$acc
  writeXStringSet(fasta, fasta_path)  
  
  junk <- tempfile()
  # setup paths for hmms and table outputs; build hmmsearch calls
  hmmsearches <- 
    tibble(hmm_path = Sys.glob(glue('{out_path}hmm/', fold, '/*'))) |> 
    mutate(
      out_path = hmm_path |> 
        str_replace('/hmm/', '/hmmsearch/') |> 
        str_replace('\\.hmm', glue('.{tag}.tbl')),
      calls = glue('hmmsearch --noali -o {junk} --tblout {out_path} {hmm_path} {fasta_path}')
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



### Models  -------------------------------------------------------------------
# update model specifications using parameters grid
make_models <- function(model, grid, name){
  models <- map(
    seq_len(nrow(grid)), ~{
      mods <- model |> 
        update(grid[.x, ])
    })
  tibble(spec = models) |> 
    mutate(model_type = name, 
           model_id = row_number()) |> 
    bind_cols(grid) |> 
    nest(params = names(grid))
}

# decision tree models
tree_models <- 
  decision_tree(mode = 'classification') |> 
  set_engine('rpart') |> 
  make_models(
    name = 'decision_tree_rpart',
    grid = grid_max_entropy(
      size = 20, tree_depth(), cost_complexity(), min_n(range = c(2L, 40L))
    ))

# logistic reg models
log_reg_models <- 
  multinom_reg(mode = 'classification') |> 
  set_engine('glmnet') |> 
  make_models(name = 'multinom_reg_glmnet',
              grid = grid_max_entropy(mixture(), penalty(), size = 20))
knn_models <- 
  nearest_neighbor(mode = 'classification') |> 
  set_engine('kknn') |> 
  make_models(name = 'nearest_neighbor',
              grid = grid_regular(neighbors(), levels = 20)
              )
# combine all model specs
models <- bind_rows(tree_models, log_reg_models, knn_models)

rm(tree_models, log_reg_models)

### CV /  Fit / Evaluate functions --------------------------------------------

#' inner_prep: function to do an inner cv fold training
#' calls prep_domains_df, build_alignments_library, build_hmm_library,
#' hmmsearch_scores, and join_hmmsearches
prep_data <- function(split, id){
  
  cat(green("\n\n------ Preparing fold: {id} --------- \n"))
  print(split)
  
  system(glue('mkdir {out_path}/align/', id))
  system(glue('mkdir {out_path}/hmm/', id))
  system(glue('mkdir {out_path}/hmmsearch/', id))
  
  # check stratification is adequate / no empty classes
  cat(green("Train and test split contain > 10 obs per subfamily?"))
  all(
    analysis(split) |> count(subfamily) |> pull(n) %>% all(.data > 10), 
    assessment(split) |> count(subfamily) |> pull(n) %>% all(.data > 10)
  ) |> 
    print()
  
  # align
  cat(green('Aligning training domains...'))
  tic()
  aligns <- 
    analysis(split) |> 
    prep_domains_df() |> 
    build_alignments_library(fold = id)
  toc()
  
  # train hmm
  cat(green('Building HMMs...'))
  hmm_check <- build_hmm_library(id)
  message('hmmbuild calls all finished without error?')
  print(all(hmm_check$hmmbuild == 0))
  
  # hmm scores test and train
  cat(green('\nScoring sequences against HMM library...'))
  tic()
  train_search <- 
    analysis(split) |> 
    hmmsearch_scores(fold = id, tag = 'train')
  test_search <- 
    assessment(split) |> 
    hmmsearch_scores(fold = id, tag = 'test')
  toc()
  
  cat(green('HMM searches all finished without error?'))
  print(all(train_search$hmmsearches == 0))
  
  # parse hmmscores and add to data frame
  train <-  
    analysis(split) |> 
    join_hmmsearches(files = train_search$out_path)
  test <-  
    assessment(split) |> 
    join_hmmsearches(files = test_search$out_path)
  
  # return tibble row with train and test data for classifiers
  cat(green('Returning datasets'))
  tibble(train = list(train), test = list(test))
}

## Fit + evaluate model on single fold
eval_model <- function(model, train, test){
  # specify evaluation functions
  class_metrics <- 
    metric_set(mcc, kap, sens, spec, precision, recall, 
               f_meas, ppv, npv, accuracy, bal_accuracy)
  
  # specify formula and scaling
  recip <- recipe(subfamily ~ ., data = train) |> 
    update_role(c('acc','description', contains('seq')), new_role = "id") |> 
    step_scale(all_predictors())
  # fit model
  fit_wf <- workflow() |>
    add_model(model) |> 
    add_recipe(recip) |>
    fit(data = train)
  # predict hold out set
  predictions <- fit_wf |> 
    predict(test |> select(-subfamily)) %>% 
    bind_cols(test)
  # collect metrics
  predictions |> 
    class_metrics(subfamily, estimate = .pred_class)
}

## TODO in progress...

## Evaluate set of models - parallelizes evaluate_model()
# 'models' needs to be tibble with column called spec with the parsnip object
# train and test are prepared tibbles from prep_data

eval_model_set <- function(models, train, test){

  plan(multisession)
  # fit / eval all models
  models |>
    mutate(metrics = future_map(
      .x = spec, .f = eval_model, train, test, 
      .options = furrr_options(packages = c("parsnip", "glmnet"), seed=NULL))
    ) |>
    unnest(metrics)
}
  
# wrapper for inner CV workflow
do_inner_cv <- function(resamples, splits, id){
  cat(magenta("\n\n***** Starting inner cv for fold  ******"))
  
  # prepare scored datasets for each fold
  prep <- resamples |> 
    mutate(prep = map2({{splits}}, {{id}}, ~prep_data(.x, .y))) |> 
    unnest(prep)
  cat(magenta("\n Evaluating Models -----"))
  
  # evaluate model set on each fold
  fold_res <- prep |> 
    mutate(rs = map2(train, test, ~eval_model_set(models, .x, .y))) |> 
    unnest(rs)
  
  # TODO eval rule-based classification
  cat(magenta("\n Summarizing results -----"))
  
  # summarize metrics across folds
  cv_res <- fold_res |> 
    group_by(model_type, model_id, .metric) |> 
    summarize(
      mean = mean(.estimate, na.rm = T),
      err = sd(.estimate, na.rm = T), 
      n_folds = sum(!is.na(.estimate)),
      values = list(.estimate),
      .groups = 'drop'
    )
  cv_res
}

# select best model hyperparameters based on metrics from inner cv
select_best_models <- function(res_summary, metric){
  cat(blue('selecting best models by {metric}'))
  res_summary |> 
    group_by(model_type) |> 
    filter(.metric == 'mcc') |> 
    arrange(desc(mean)) |> 
    slice_head(n = 1) |> 
    select(model_type, model_id) |> 
    left_join(models, by = c("model_type", "model_id")) |> 
    left_join(res_summary, by = c("model_type", "model_id")
    ) |> 
    ungroup()
}


# wrapper for outer cv
# TODO eval rule-based classification

fit_nested_cv <- function(nestcv){
  
  cat(yellow$bold('______      Beginning inner CVs       ________'))
  # do inner cv for each outer fold
  results_df <- nest_cv |> 
    mutate(
      inner_cv = map(inner_resamples, 
                     ~do_inner_cv(.x, .x$inner_splits, .x$inner_id)),
    ) 
  
  cat(yellow$bold('______     Selecting model tunings    ________'))
  # select best model tunings for each outer fold
  results_df <- results_df |> 
    mutate(best_tune = map(inner_cv, ~{
      select_best_models(.x, metric = 'mcc') |> 
        select(model_type, model_id, spec) |>
        distinct()}
    ))  
  
  cat(yellow$bold('______   Preparing outer split data   ________'))
  # prepare outer splits
  results_df <- results_df |> 
    mutate(prep = map2(outer_splits, outer_id, ~prep_data(.x, .y))) |> 
    unnest(prep)  
  
  cat(yellow$bold('______    Evaluating tuned models      ________'))
  # evaluate tuned models on outer splits
  results_df <- results_df |>
    mutate(results = pmap(
      .l = list(best_tune, train, test),
      .f = function(best_tune, train, test){
        eval_model_set(models = best_tune,  train = train, test = test)
      })
    )
  cat(yellow$bold('______           Finished             ________'))
  return(results_df)
}

# MAIN / TRAIN & TEST ----------------------------------------

## TODO max & threshold classifier
## TODO full dataset, minimal resampling...


nest_cv_results <- nest_cv |> fit_nested_cv()

beepr::beep()

results_df |> 
  select(results) |> 
  unnest(cols = c(results)) |> 
  group_by(model_type, model_id, .metric) |> 
  summarize(
    mean = mean(.estimate, na.rm = T),
    err = sd(.estimate, na.rm = T), 
    n_folds = sum(!is.na(.estimate)),
    values = list(.estimate),
    .groups = 'drop'
  ) |> 
  filter(.metric %in% c('mcc')) |> 
  unnest(values)


# -----
# 

# 
# # do inner cv for each outer fold
# inner_cvs <-
#   nest_cv |> 
#   pull(inner_resamples) |> 
#   map(~do_inner_cv(.x, splits = .x$inner_splits,  id = .x$inner_id)) 
# 
# # select best model tunings for each outer fold
# nest_cv_mods <- inner_cvs |> 
#   map(~select_best_models(.x, metric = 'mcc')) |> 
#   map(~{.x |> 
#       select(model_type, model_id, spec) |>
#       distinct()
#   }) 
# beepr::beep()
# 
# 
# # attach to each fold: the best tuning specs in inner CV for each model 
# nest_cv_tuned <- nest_cv |> mutate(mods = nest_cv_mods)
# 
# # prep all outer splits
# tic()
# prep <- nest_cv_tuned |> 
#   mutate(prep = map2(outer_splits, outer_id, ~prep_data(.x, .y))) |> 
#   unnest(prep)
# toc()
# prep
# 
# # evaluate models on outer resamples
# tic()
# outer_cv_res <- prep %>%
#   mutate(rs = pmap(
#     .l = list(mods, train, test), 
#     .f = function(mods, train, test){
#       eval_model_set(models = mods,  train = train, test = test)
#     })
#   ) 
# toc()
# 
# return(outer_cv_res)



# results_df <- 
#   nest_cv |> 
#   # do inner cv for each outer fold
#   mutate(
#     inner_cv = map(inner_resamples, 
#                    ~do_inner_cv(.x, .x$inner_splits, .x$inner_id)),
#     )  |> 
#   # select best model tunings for each outer fold
#   mutate(best_tune = map(inner_cv, ~{
#     select_best_models(.x, metric = 'mcc') |> 
#       select(model_type, model_id, spec) |>
#       distinct()}
#   ))  |> 
#   mutate(prep = map2(outer_splits, outer_id, ~prep_data(.x, .y))) |> 
#   unnest(prep)  |>
#   mutate(results = pmap(
#     .l = list(best_tune, train, test),
#     .f = function(best_tune, train, test){
#       eval_model_set(models = best_tune,  train = train, test = test)
#     })
#   )