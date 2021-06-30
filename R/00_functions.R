### Functions for nested_cv

### Libraries ---------------------------------------------------------------
library(tidyverse)
library(tidymodels)
library(Biostrings)
library(furrr)
library(DECIPHER)
library(glue)
library(here)
library(tictoc)
library(crayon)

### Prep Functions ---------------------------------------

## Prepare domains for alignment/hmm 
prep_domains_df <- function(training){
  training |> 
    # subset domain subfamilies
    filter(subfamily !=  'non_integrase') |> 
    select(subfamily, acc, dom_seq) |> 
    # remove any duplicated domain sequences
    group_by(dom_seq) |> sample_n(1) |> ungroup() |> 
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
  plan(multisession, workers = availableCores())
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
# takes protein seq column and makes temporary fasta file to send to hmmer
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


#' prep_data:
#' wrapper function to do alignments, hmmbuilds, hmmsearches, & parse searches
#' does one train / test split 
#' calls prep_domains_df, build_alignments_library, build_hmm_library,
#' hmmsearch_scores, and join_hmmsearches
prep_data <- function(split, id){
  
  cat('\n', green$bold(glue("Preparing fold: {id}")))
  # create directories
  system(glue('mkdir {out_path}/align/', id))
  system(glue('mkdir {out_path}/hmm/', id))
  system(glue('mkdir {out_path}/hmmsearch/', id))
  
  # check stratification is adequate / no empty classes
  strat_test <- all(
    analysis(split) |> count(subfamily) |> pull(n) %>% all(.data > 10), 
    assessment(split) |> count(subfamily) |> pull(n) %>% all(.data > 10)
  )
  cat('\n', white(glue("Train and test split contain > 10 obs per subfamily? {strat_test}")))
  
  # align
  cat('\n', white('Aligning training domains...'))
  tic()
  aligns <- 
    analysis(split) |> 
    prep_domains_df() |> 
    build_alignments_library(fold = id)
  toc()
  
  # train hmm
  cat(white(' Building HMMs...'))
  hmm_check <- build_hmm_library(id)
  hmmbuild_test <- all(hmm_check$hmmbuild == 0)
  cat('\n', white(glue('hmmbuild calls all finished without error? {hmmbuild_test}')))
  
  # hmm scores test and train
  cat('\n', white('Scoring sequences against HMM library...'))
  tic()
  train_search <- 
    analysis(split) |> 
    hmmsearch_scores(fold = id, tag = 'train')
  test_search <- 
    assessment(split) |> 
    hmmsearch_scores(fold = id, tag = 'test')
  toc()
  
  hmmsearch_test <- all(train_search$hmmsearches == 0)
  cat(white(glue(' HMM searches all finished without error? {hmmsearch_test}')))
  
  # parse hmmscores and add to data frame
  train <-  
    analysis(split) |> 
    join_hmmsearches(files = train_search$out_path)
  test <-  
    assessment(split) |> 
    join_hmmsearches(files = test_search$out_path)
  
  # return tibble row with train and test data for classifiers
  cat(white(' Returning datasets'))
  tibble(train = list(train), test = list(test))
}

### Threshold classifier --------------------------------------------------

set_thresholds <- function(train_scores){
  train_scores |> 
    select(subfamily, Arch1:Xer) |> 
    pivot_longer(cols = Arch1:Xer,  names_to = 'hmm_name', 
                 values_to = 'threshold') |> 
    filter(subfamily == hmm_name) |> 
    group_by(subfamily) |> 
    filter(threshold == min(threshold)) |> 
    ungroup() |> 
    distinct() |> 
    select(-subfamily)
}

hmm_threshold_class <- function(test, thresholds){
  # for pred_class factor to have the same levels
  levs <- c("Arch1",  "Arch2", "Candidate", "Cyan", "IntKX", "Int_BPP-1",
            "Int_Brujita", "Int_CTnDOT", "Int_Des", "Int_P2", "Int_SXT", 
            "Int_Tn916" ,"Integron", "Myc", "RitA", "RitB", "RitC", "TnpA", 
            "TnpR", "Xer", "non_integrase")
  test |> 
    select(acc, Arch1:Xer) |> 
    pivot_longer(cols = Arch1:Xer, 
                 names_to = 'hmm_name', 
                 values_to = 'hmm_score') |> 
    group_by(acc) |> 
    filter(hmm_score == max(hmm_score)) |> 
    ungroup() |> 
    left_join(thresholds, by = "hmm_name") |> 
    mutate(hmm_name = ifelse(hmm_score < threshold, 'non_integrase', hmm_name)) |>     select(-threshold) |> 
    distinct() |> 
    ungroup() |> 
    transmute(acc, .pred_class = factor(hmm_name, levels = levs)) |> 
    left_join(test, by = 'acc')
}

# sets thresholds from training, predict test, evaluate metrics
eval_threshold_classififer <- function(train, test){
  # set thresholds
  thresholds <- set_thresholds(train)
  # predict based on top score & threshold
  predictions <- hmm_threshold_class(test, thresholds)
  # collect metrics
  class_metrics <- metric_set(mcc, kap, sens, spec, precision, recall, 
                              f_meas, ppv, npv, accuracy, bal_accuracy)
  predictions |>
    class_metrics(subfamily, estimate = .pred_class) 
}

### CV /  Fit / Evaluate functions --------------------------------------------

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
    predict(object = test |> select(-subfamily)) %>% 
    bind_cols(test)
  # collect metrics
  predictions |> 
    class_metrics(subfamily, estimate = .pred_class)
}



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
  
  cat("\n", blue$bold("Starting inner CV"), '\n')
  # prepare scored datasets for each fold
  prep <- resamples |> 
    mutate(prep = map2({{splits}}, {{id}}, ~prep_data(.x, .y))) |> 
    unnest(prep)
  
  cat('\n', blue$bold("Evaluating models"))
  tic()
  # evaluate model set on each fold
  fold_res <- prep |> 
    mutate(rs = map2(train, test, ~eval_model_set(models, .x, .y))) |> 
    unnest(rs)
  toc()
  
  cat('\n', blue$bold("Summarizing results"))
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
  cat('\n', blue('selecting best models by {metric}'), '\n')
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
fit_nested_cv <- function(nestcv){
  
  cat(yellow$bold$underline('Beginning inner CVs'))
  # do inner cv for each outer fold
  results_df <- nest_cv |> 
    mutate(
      inner_cv = map(inner_resamples, 
                     ~do_inner_cv(.x, .x$inner_splits, .x$inner_id)),
    ) 
  
  cat('\n', yellow$bold$underline('Selecting model tunings'))
  # select best model tunings for each outer fold
  results_df <- results_df |> 
    mutate(best_tune = map(inner_cv, ~{
      select_best_models(.x, metric = 'mcc') |> 
        select(model_type, model_id, spec) |>
        distinct()}
    ))  
  
  cat('\n', yellow$bold('Preparing outer split data'))
  # prepare outer splits
  results_df <- results_df |> 
    mutate(prep = map2(outer_splits, outer_id, ~prep_data(.x, .y))) |> 
    unnest(prep)  
  
  cat('\n', yellow$bold('Evaluating tuned models'))
  # evaluate tuned models on outer splits
  results_df <- results_df |>
    mutate(results = pmap(
      .l = list(best_tune, train, test),
      .f = function(best_tune, train, test)
        eval_model_set(models = best_tune,  train = train, test = test)
    )
    )
  cat('\n', yellow$bold('Finished'))
  return(results_df)
}



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


