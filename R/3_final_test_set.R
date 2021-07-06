
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
source('./R/00_get_model_specs.R')

set.seed(12345)

folder <- 'final_model_07-04'
system(glue('mkdir ./results/{folder}/'))
system(glue('mkdir ./results/{folder}/align/'))
system(glue('mkdir ./results/{folder}/hmm/'))
system(glue('mkdir ./results/{folder}/hmmsearch/'))

train <- read_rds('./data/classif_train_set.rds')
test <- read_rds('./data/classif_test_set.rds')


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


## Main: Prep ----
tic()
# do alignments
cat(white(' Doing alignments ... \n\n'))
aligns <- 
  prep_domains_df(train) |> 
  build_alignments_library2(save_path = glue('./results/{folder}/align/'))
toc()


# train hmms
cat(white(' Building HMMs ...\n\n'))
hmm_check <- build_hmm_library2(msa_path =  glue('./results/{folder}/align/'),
                                hmm_path = glue('./results/{folder}/hmm/'))

hmmbuild_test <- all(hmm_check$hmmbuild == 0)
cat('\n', white(glue('hmmbuild calls all finished without error?\n {hmmbuild_test}')))


# generate hmm scores for test and train (saves to outfiles)
cat('\n', white(' Scoring sequences against HMM library...'))
tic()
train_scores <-
  hmmsearch_scores2(train, hmm_path = glue('./results/{folder}/hmm/'), tag = 'train')
test_scores <-
  hmmsearch_scores2(test, hmm_path = glue('./results/{folder}/hmm/'), tag = 'test')
toc()

# parse hmmscores and add to data frame
train_prep <- train |> join_hmmsearches(files = train_scores$out_path)
test_prep <- test |> join_hmmsearches(files = test_scores$out_path)

## prepared data
prep <- tibble(train_prep = list(train_prep), 
                test_prep = list(test_prep))

write_rds(prep, './results/final_prep.rds')
write_rds(aligns, './results/final_train_aligns.rds')

rm(train, test, train_scores, test_scores, hmm_check)


prep <- read_rds('./results/final_prep.rds')

train_prep <- prep$train[[1]]
test_prep <- prep$test[[1]]


# specify modelling workflow with scaling and ignore seqs and id
recip <- recipe(subfamily ~ ., data = train_prep) |> 
  update_role(c('acc', contains('seq')), new_role = "id") |> 
  step_scale(all_predictors())

# see scaled data
recip |> prep() |> juice()

wkfl <- workflow() |> 
  add_recipe(recip)
wkfl


# create a workflow for model spec with recipe, then fit to training set
fit_workflow <- function(spec, wkfl, train_prep) {
 return(wkfl |>  add_model(spec) |>  fit(data = train_prep))
}

# apply fitting to all models
plan(multisession, workers = 8)

fitted_models <- models |> 
  mutate(fitted_wkfl = future_map(
    .x = spec, 
    .f = ~fit_workflow(.x, wkfl, train_prep),
    .options = furrr_options(
      packages = c("parsnip","tidymodels", "glmnet","rpart", "kknn"),
      seed = NULL)
    ))
beepr::beep()

fitted_models
fitted_models$fitted_wkfl[[1]]
fitted_models$fitted_wkfl[[21]]
fitted_models$fitted_wkfl[[41]]
fitted_models$fitted_wkfl[[43]]


# predict test data
get_predictions <- function(fitted_wkfl, test){
  predictions <- 
    fitted_wkfl |>
    predict(test |> select(-subfamily)) %>%
    bind_cols(test |> select(subfamily))
}

# apply get_predictions to each fitted workflow
plan(multisession, workers = 8)

fitted_models_w_preds <-  fitted_models |> 
  mutate(preds = future_map(
    .x = fitted_wkfl, 
    .f = ~get_predictions(.x, test_prep),
    .options = furrr_options(
      packages = c("parsnip","tidymodels", "glmnet","rpart", "kknn"),
      seed = NULL)
    ))

beepr::beep()

# specify evaluation functions
class_metrics <- 
  metric_set(mcc, kap, sensitivity, specificity, precision, recall, 
             f_meas, ppv, npv, accuracy, bal_accuracy)




## Evaluate prediction performance



## Fit/evaluate rule-based classifier
thresh_res <- prep |>
  mutate(threshold_res = map2(train, test, ~eval_threshold_classififer(.x, .y))) |>
  select(threshold_res) |>
  unnest(threshold_res) |>
  mutate(model_type = 'score & threshold rule',
         model_id = NA)

# final_test_res <- 
#   bind_rows(thresh_res, metrics) |> 
#   mutate(model_id = model_id)
# rm(thresh_res)
# 
# final_test_res |> View()
# 
# 
# 
# 
# 
# best_in_cv <- 
#   read_rds('./results/3x3_regular_CV_07-02/best_models.rds')
# 
# best_in_cv |> 
#   left_join(final_test_res)
#   
# 
# 
# final_test_res |> 
#   ggplot(aes(x = estimate, y = fct_rev(model_type)))
