
library(tidyverse)
library(tidymodels)
library(Biostrings)
library(furrr)
library(DECIPHER)
library(glue)
library(here)
library(tictoc)
library(crayon)

source('./R/00_functions.R')
models <- read_rds('data/unfitted_parsnip_model_set.rds')
options(tibble.print_max = 50)

set.seed(12345)

folder <- 'final_model_07-11'
system(glue('mkdir ./results/{folder}/'))
system(glue('mkdir ./results/{folder}/align/'))
system(glue('mkdir ./results/{folder}/hmm/'))
system(glue('mkdir ./results/{folder}/hmmsearch/'))
outpath <- glue('./results/{folder}')
# split_id <- 'full_training'

train <- read_rds('./data/classif_train_set.rds')
test <- read_rds('./data/classif_test_set.rds')

best_cv_models <- 
  read_rds('./results/regular_cv_07-10/best_models.rds')

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
  wkfl |>  
    add_model(spec) |>  
    fit(data = train_prep)
}

# predict test data
get_predictions <- function(fitted_wkfl, test){
  fitted_wkfl |>
    predict(test |> select(-subfamily))
}

# specify evaluation functions
collect_pref_metrics <- function(preds, truth, met_set){
  preds <- preds$.pred_class
  my_metrics <- metric_set(mcc, kap, sensitivity, 
                           specificity, precision, recall, 
                           f_meas, ppv, npv, accuracy, bal_accuracy)
  tibble(t=truth, p=preds) |> 
    my_metrics(estimate = p, truth = t)
}


# Main ----

### Prep ----


# do alignments
tic()
cat(white(' Doing alignments ... \n\n'))
aligns <- 
  prep_domains_df(train) |> 
  build_alignments_library2(save_path = glue('./results/{folder}/align/'))
toc()
beepr::beep()

# train hmms
cat(white(' Building HMMs ...\n\n'))
hmm_check <- build_hmm_library2(
  msa_path =  glue('./results/{folder}/align/'),
  hmm_path = glue('./results/{folder}/hmm/')
  )
hmm_check
hmmbuild_test <- all(hmm_check$hmmbuild == 0)
cat('\n', white(glue('hmmbuild calls all finished without error?\n {hmmbuild_test}')))
beepr::beep()

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
beepr::beep()

# parse hmmscores and add to data frame
train_prep <- train |> join_hmmsearches(files = train_scores$out_path)
test_prep <- test |> join_hmmsearches(files = test_scores$out_path)

## prepared data
prepped_data <- tibble(train_prep = list(train_prep), 
                       test_prep = list(test_prep))
prepped_data

write_rds(prepped_data, glue('./results/{folder}/prepped_data.rds'))
rm(train, test, train_scores, test_scores, hmm_check)


### Modelling --------------------------------------------------------

prepped_data <- read_rds(glue('./results/{folder}/prepped_data.rds'))

train_prep <- prepped_data$train_prep[[1]] |> 
  select(-contains('seq'))

test_prep <- prepped_data$test_prep[[1]] |> 
  select(-contains('seq'))

rm(prepped_data)

train_prep |> count(subfamily)
test_prep |> count(subfamily)

# TODO change recipe
# recip <- 
#   recipe(subfamily ~ ., data = train_prep) |> 
#   update_role(c('acc', contains('seq')), new_role = "id") |> 
#   step_(all_predictors())

# specify modelling workflow with scaling and ignore seqs and id
recip <- 
  recipe(subfamily ~ ., data = train) |> 
  update_role(acc, new_role = "id") |> 
  step_smote(subfamily, over_ratio = 0.25) |> 
  step_normalize(all_predictors())

# see scaled data
recip |> prep() |> juice() |> skimr::skim()

wkfl <- workflow() |> 
  add_recipe(recip)
wkfl

# selected models from 2c_regular_CV.R
prev_best_models <- 
  bestbest_cv_models |> 
  filter(!is.na(model_id))  |>
  nest(reg_cv_res = c(.metric, mean, err, n_folds, values)) |> 
  left_join(models, by = c("model_type", "model_id"))


#### Fit, predict, evaluate -----

cat('\n', white('Fitting, predicting, and evaluating performance'))
plan(multisession, workers = 8)
fitted_models <- 
  models |> 
  mutate(
    # fit each model to the prepared training data
    fitted_wkfl = future_map(
      .x = spec, 
      .f = ~fit_workflow(.x, wkfl, train_prep),
      .options = furrr_options(
        packages = c("parsnip","tidymodels", "glmnet","rpart", "kknn"),
        seed = NULL)
    ),
    # get predictions from each model for final test set
    preds = future_map(fitted_wkfl, ~get_predictions(.x, test_prep)),
    # evaluate prediction performance
    final_metrics = future_map(preds, ~collect_pref_metrics(preds = .x, truth = test_prep$subfamily))
    ) |> 
  select(-spec)

beepr::beep()


write_rds(fitted_models, glue('./results/{folder}/fitted_models.rds'))

best_models |> 
  left_join(fitted_models |> select(model_type, model_id, final_metrics))

# best_models |> 
#   unnest(reg_cv_res)
# models |>
#   filter(model_type == 'multinom_reg_glmnet') |>
#   unnest(params) |>
#   ggplot(aes(penalty, mixture)) +
#   geom_text(aes(label = model_id)) +
#   scale_x_log10()

fitted_models |> 
  select(model_type, model_id, final_metrics) |> 
  unnest(final_metrics) |> 
  filter(.metric == 'mcc') |> 
  ggplot(aes(x = .estimate, y = fct_rev(factor(model_id)), fill = model_type)) +
  facet_wrap(~model_type) +
  geom_col()
  

# wtf, now model 7 is terrible?
fitted_models |> 
  filter(model_type == 'multinom_reg_glmnet' & model_id == 7) |> 
  select(preds) |> unnest(cols = c(preds)) |> 
  bind_cols(truth = test_prep$subfamily) |> 
  my_metrics(truth = truth, estimate = .pred_class, na_rm = T)
  


### Eval rule-based classifier ----

prepped_data <- read_rds(glue('./results/{folder}/final_prep.rds'))

thresh_res <- prepped_data |>
  mutate(threshold_res = map2(train_prep, test_prep, ~eval_threshold_classififer(.x, .y))) |>
  select(threshold_res) |>
  unnest(threshold_res) |>
  mutate(model_type = 'score & threshold rule', model_id = NA)

thresh_res

final_validation_test_results <- 
  fitted_models |> 
  select(model_type, model_id, final_metrics) |> 
  bind_rows(thresh_res |> nest(final_metrics = c('.metric', '.estimator', '.estimate')))

write_rds(final_validation_test_results,
          './results/{folder}/final_validation_results.rds')

