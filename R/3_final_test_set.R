
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


# hmm scores test and train
cat('\n', white(' Scoring sequences against HMM library...'))
tic()
train_prep <-
  hmmsearch_scores2(train, hmm_path = glue('./results/{folder}/hmm/'), tag = 'train')
test_prep <-
  hmmsearch_scores2(test, hmm_path = glue('./results/{folder}/hmm/'), tag = 'test')
toc()

# parse hmmscores and add to data frame
train <- train |> join_hmmsearches(files = train_prep$out_path)
test <- test |> join_hmmsearches(files = test_prep$out_path)

## prepared data
prep <- tibble(train = list(train), test = list(test))
write_rds(prep, './final_prep.rds')

prep

## Fit/evaluate models
res <- prep |> 
  mutate(rs = map2(train, test, ~ eval_model_set(models, .x, .y)))


