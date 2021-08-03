## 3 Final validation on holdout set

library(tidyverse)
library(tidymodels)
library(Biostrings)
library(furrr)
library(DECIPHER)
library(glue)
library(here)
library(tictoc)
library(crayon)
library(themis)
library(beepr)
options(tibble.print_max = 50)

# data prep and modelling functions
source('./R/00_functions.R')

# where in ./results/ to save files to?
folder <- 'final_model_07-11'

# make directories
system(glue('mkdir ./results/{folder}/'))
system(glue('mkdir ./results/{folder}/align/'))
system(glue('mkdir ./results/{folder}/hmm/'))
system(glue('mkdir ./results/{folder}/hmmsearch/'))
outpath <- glue('./results/{folder}')

## Functions ----

# Wrapper to apply align_domains to each subfamily
build_alignments_library2 <- function(domains, save_path){
  # do align_domains in parallel
  plan(multisession, workers = availableCores())
  domains |> mutate(
    aligned = future_map(
      .x = data, 
      .f = ~ align_domains(.x, dest = save_path),
      .options = furrr_options(scheduling = Inf)
    )
  )
}

# create paths, build call, do calls to hmmbuild for each subfamily
build_hmm_library2 <- function(msa_path, hmm_path){
  
  temp <- tempfile()
  plan(multisession, workers = availableCores())
  
  tibble(msa_path = Sys.glob(glue('{msa_path}/*'))) |> 
    mutate(
      hmm_path = str_replace_all(msa_path, 'align|aln', 'hmm'),
      hmmbuild_call = glue('hmmbuild -o {temp} {hmm_path} {msa_path}'),
      hmmbuild = future_map_dbl(hmmbuild_call, ~ system(.x))
    )
}

# score a dataframe of sequences against the HMM library with hmmsearch
# takes protein seq column and makes temporary fasta file to send to hmmer
hmmsearch_scores2 <- function(df, hmm_path, tag){
  
  # make temp fasta file of seqs to score against hmm library
  fasta_path <- tempfile()
  fasta <- Biostrings::AAStringSet(df$prot_seq)
  names(fasta) <- df$acc
  writeXStringSet(fasta, fasta_path)  
  
  junk <- tempfile()
  # setup paths for hmms and table outputs; build hmmsearch calls
  hmmsearches <- 
    tibble(hmm_path = Sys.glob(glue('{hmm_path}*'))) |> 
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


# create a workflow for model spec with recipe, then fit to training set
fit_workflow <- function(spec, wkfl, train_prep) {
  wkfl |> add_model(spec) |> fit(data = train_prep)
}

# predict test data
get_predictions <- function(fitted_wkfl, test){
  fitted_wkfl |> predict(test |> select(-subfamily))
}

# specify evaluation functions
collect_perf_metrics <- function(preds, truth, met_set){
  preds <- preds$.pred_class
  my_metrics <- metric_set(mcc, kap, sensitivity, 
                           specificity, precision, recall, 
                           f_meas, ppv, npv, accuracy, bal_accuracy)
  tibble(t=truth, p=preds) |> 
    my_metrics(estimate = p, truth = t)
}


# Main ------------------------------------------------------------

### Data ----

train <- read_rds('./data/classif_train_set.rds')
train
test <- read_rds('./data/classif_test_set.rds')
test

# unfitted model specifications
models <- 
  read_rds('data/unfitted_parsnip_model_set.rds') |> 
  filter(model_type != 'decision tree')
models

# best models from CV
selected_models <-
  read_rds('./results/regular_cv_07-10/best_models.rds') |> 
  filter(!is.na(model_id))  |>
  nest(reg_cv_res = c(.metric, mean, err, n_folds, values)) |>
  left_join(models, by = c("model_type", "model_id"))
selected_models
rm(models)

### Prep --------------------------------------------------------

# do alignments
tic()
cat(white(' Doing alignments ... \n\n'))
aligns <- 
  prep_domains_df(train) |> 
  build_alignments_library2(save_path = glue('./results/{folder}/align/'))
toc()
beep()

# train hmms
cat(white(' Building HMMs ...\n\n'))
hmm_check <- build_hmm_library2(
  msa_path =  glue('./results/{folder}/align/'),
  hmm_path = glue('./results/{folder}/hmm/')
  )
hmm_check
hmmbuild_test <- all(hmm_check$hmmbuild == 0)
cat('\n', white(glue('hmmbuild calls all finished without error?\n {hmmbuild_test}')))
beep()

# generate hmm scores for test and train (saves to outfiles)
cat('\n', white(' Scoring sequences against HMM library...'))
tic()
train_scores <-
  hmmsearch_scores2(train, 
                    hmm_path = glue('./results/{folder}/hmm/'),
                    tag = 'train')
test_scores <-
  hmmsearch_scores2(test, 
                    hmm_path = glue('./results/{folder}/hmm/'),
                    tag = 'test')
toc()
beep()

# parse hmmscores and add to data frame
train_prep <- train |> join_hmmsearches(files = train_scores$out_path)
test_prep <- test |> join_hmmsearches(files = test_scores$out_path)

## prepared data - scored but not normalized or SMOTE'd
prepped_data <- tibble(train_prep = list(train_prep), 
                       test_prep = list(test_prep))
prepped_data

## save prep work
write_rds(prepped_data, glue('./results/{folder}/prepped_data.rds'))
rm(train, test, train_scores, test_scores, hmm_check)


### Modelling ------------------------------------------


# take training data from prep
train_prep <-
  read_rds(glue('./results/{folder}/prepped_data.rds')) |> 
  pull(train_prep) |> 
  pluck(1) |> 
  select(-contains('seq'))

train_prep
train_prep |> count(subfamily)

#### Workflows -------------------------------------

# specify modelling workflow with scaling and ignore seqs and id
recip <- 
  recipe(subfamily ~ ., data = train_prep) |> 
  update_role(acc, new_role = "id") |> 
  step_smote(subfamily, over_ratio = 0.25) |> 
  step_normalize(all_predictors())

# see scaled data
recip |> prep() |> juice() |> skimr::skim()

# create workflow for modelling
wkfl <- workflow() |> add_recipe(recip)
wkfl


#### Fit & Eval ------------------------------------------------

# Fit best CV models and predict validation set and evaluate metrics

# load prepared test split data
test_prep <- 
  read_rds(glue('./results/{folder}/prepped_data.rds')) |> 
  pull(test_prep) |> 
  pluck(1) |> 
  select(-contains('seq'))

# fit each model to the prepared training data
tic()
plan(multisession, workers = 8)

fit_models <- selected_models |>
  mutate(
    fitted_wkfl = future_map(
      .x = spec,
      .f = ~fit_workflow(.x, wkfl, train_prep),
      .options = furrr_options(
        packages = c('parsnip','tidymodels','ranger', 'kknn',
                     'glmnet', 'themis'),
        seed = TRUE)
    ),
    preds = map(
      .x = fitted_wkfl, 
      .f = ~get_predictions(fitted_wkfl = .x, test = test_prep)
    ),
    # evaluate prediction performance
    final_metrics = map(
      .x = preds, 
      .f = ~collect_perf_metrics(preds = .x, 
                                 truth = test_prep$subfamily)
    )) |> 
  select(-spec)
  
toc()
beep()

write_rds(x = fit_models,
          file = glue('./results/{folder}/fit_models.rds'),
          compress = 'gz')


### Eval rule-based classifier ----

prepped_data <- read_rds(glue('./results/{folder}/prepped_data.rds'))

thresh_res <- 
  prepped_data |>
  mutate(threshold_res = map2(
    .x = train_prep, 
    .y = test_prep, 
    .f = ~eval_threshold_classififer(.x, .y)
  )) |>
  select(threshold_res) |>
  unnest(threshold_res) |>
  mutate(model_type = 'score & threshold rule', model_id = NA)

thresh_res
rm(prepped_data)

## old_TODO juice recip and test classifier *
## taken care of in 4b_plots

normalized_smoted <- recip
normalized_threshold_res <- 
  recip |> prep() |> juice() |> 
  mutate(threshold_res = map2(
    .x = train_prep, 
    .y = test_prep, 
    .f = ~eval_threshold_classififer(.x, .y)
  )) |>
  select(threshold_res) |>
  unnest(threshold_res) |>
  mutate(model_type = 'score & threshold rule', model_id = NA)

# 
# thresh_res
# rm(prepped_data)

## Save Final Results  ----

final_res <- fit_models |>
  select(model_type, model_id, final_metrics) |>
  bind_rows(
    thresh_res |>
      nest(final_metrics = c('.metric', '.estimator', '.estimate'))
  )

final_res |> unnest(final_metrics) |> print(n=100)

write_rds(
  x = final_res,
  file = glue('./results/{folder}/final_validation_results.rds'),
  compress = 'gz'
)


## Check performance -------


## Finished

