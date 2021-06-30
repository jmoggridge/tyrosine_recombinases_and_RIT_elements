# 3_CV

dir <- 'classifier_02/'

## Combine datasets  ------

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

set.seed(123)
# combine datasets
full_dataset <- 
  bind_rows(smart_df, non_integrases) |> 
  # TODO (remove downsampling)
  # downsample to 10000
  group_by(subfamily) |> 
  # TODO make downsampling consistent with nested cv n=15k
  slice_sample(n = 100, replace = F) |> 
  ungroup()

rm(smart_df, non_integrases)

## Initial split -------------------------------------------------

set.seed(54321)

df_split <- initial_split(full_dataset, 0.75, strata = subfamily)
train <- training(df_split)
test <- testing(df_split)

rm(full_dataset)


## 2x2-fold CV

cv_splits <- vfold_cv(train, v = 3, repeats = 2, strata = subfamily)

fit_cv <- function(cv_splits){
  # prep
  cv_splits <- cv_splits |> 
    mutate(prep = map2(outer_splits, outer_id, ~prep_data(.x, .y))) |> 
    unnest(prep)  
  
  # evaluate models
  
}


