
# need these packages
library(tidyverse)
library(tidymodels)
library(furrr)
library(glue)
library(Biostrings)

# eventually will replace this script with a package
source('R/predict_for_others-pkg_beginnings/predict_integrase_subfamilies.R')

# Maverick's data
mav <- 
  # create dataframe of filepaths to faa files
  tibble(files = Sys.glob('data/Maverick_MGEs_for_integrase_analysis/*')) |> 
  mutate(
    # read fasta files
    fasta = map(files, ~Biostrings::readAAStringSet(.x)),
    # parse names and sequences
    name = map(fasta, names),
    seq = map(fasta, paste),
    file = str_remove(files, 'data/Maverick_MGEs_for_integrase_analysis/')
    ) |> 
  select(file, name, seq) |> 
  unnest(cols = c(name, seq))
  
glimpse(mav)

# Get integrase subfamily predictions
mav_preds <- 
  predict_integrases(mav, names = name, seqs = seq) 

# no classifications without support from 2 different models
mav_preds |> 
  filter(consensus_pred == 'no consensus') |> 
  nrow()
  
# bar chart
mav_preds |> ggplot(aes(y = consensus_pred)) + geom_bar()

# keep just the integrase hits, discard other & no consensus
integrases <- 
  mav_preds |> 
  select(-c(Arch1:Xer, matches('rf|knn|glmnet'))) |> 
  filter(consensus_pred != 'Other') |> 
  transmute(file, name, int_subfamily = consensus_pred)

write_csv(integrases, '~/Desktop/Maverick_integrases.csv')
