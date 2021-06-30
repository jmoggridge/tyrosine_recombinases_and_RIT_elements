library(tidyverse)
library(furrr)
library(tidymodels)

#' Precision: proportion of correct positive identification (TP / TP + FP)
#' Recall: proportion of actual positives predicted correctly (TP / TP + FN)
#' 
#' 
#' 
#' 
#' 
#' 

# TODO finish rules-based classification

train <- read_rds('./data/train_df_scored.rds')
train
test <- read_rds('./data/test_df_scored.rds') 
test


test |>  count(subfamily) |> arrange(desc(n)) |> print.AsIs()

## No NA's
# sum(is.na(test |> select(Arch1:Xer)))
# sum(is.na(train |> select(Arch1:Xer)))
# sum(is.na(test |> select(AA:YY)))
# sum(is.na(train |> select(AA:YY)))

# set threshold for each model as the lowest score of true members
gathering_thresholds <- 
  train |> 
  select(subfamily, Arch1:Xer) |> 
  pivot_longer(cols = Arch1:Xer, names_to = 'hmm_name', values_to = 'threshold') |> 
  filter(subfamily == hmm_name) |> 
  group_by(subfamily) |> 
  filter(threshold == min(threshold)) |> 
  ungroup() |> 
  distinct() |> 
  select(-subfamily)

print.data.frame(gathering_thresholds)
  

## rules-based classification
test_preds <- 
  test |> 
  select(acc, subfamily, Arch1:Xer) |> 
  pivot_longer(cols = Arch1:Xer, names_to = 'hmm_name', values_to = 'hmm_score') |> 
  group_by(acc, subfamily) |> 
  filter(hmm_score == max(hmm_score)) |> 
  left_join(gathering_thresholds, by = "hmm_name") |> 
  mutate(hmm_name = ifelse(hmm_score < threshold, 'non_integrase', hmm_name)) |> 
  select(-threshold) |> 
  distinct() |> 
  ungroup() |> 
  ## actual and predicted classes need to have the same ordering of factor levels
  arrange(subfamily) |> 
  mutate(subfamily = as_factor(subfamily)) |> 
  arrange(hmm_name) |> 
  mutate(predicted = as_factor(hmm_name)) 


# which integrases did not meet their own model's cutoffs?
test_preds |> 
  filter(predicted == 'non_integrase') |> 
  count(subfamily, predicted)

  
test_preds |> 
  filter(predicted == 'Int_Des')



## confusion matrix plot
test_preds |> 
  count(subfamily, hmm_name) |> 
  left_join(
    test |> count(subfamily, name = 'count')
    ) |> 
  mutate(predicted = hmm_name, 
         actual = fct_rev(subfamily)) |> 
  ggplot(aes(predicted, actual)) +
  geom_point(aes(size = n)) +
  scale_size(trans = 'log10') +
  theme_light() +
  theme(axis.text.x.bottom = element_text(angle = -90,hjust = 0, vjust = 0)) +
  labs(subtitle = 'Gathering-threshold based classification:')


# conf_matrix table  
cm <- table(test_preds$subfamily, test_preds$predicted)

# precision
prec_ls <- 
  (diag(cm) / rowSums(cm)) |> 
  enframe(name = 'subfamily', value = 'precision')
# recall
recall_ls <- 
  (diag(cm) / colSums(cm)) |>
  enframe(name = 'subfamily', value = 'recall')

class_wise_pr <- full_join(prec_ls, recall_ls, by = "subfamily") 
class_wise_pr |> print.AsIs()

get_measures <- function(df, truth, estimate){
  bind_rows(
    df |> accuracy({{truth}}, {{estimate}}),
    df |> bal_accuracy({{truth}}, {{estimate}}),
    df |> kap({{truth}}, {{estimate}}),
    df |> f_meas({{truth}}, {{estimate}}),
    df |> mcc({{truth}}, {{estimate}}),
    df |> spec({{truth}}, {{estimate}}),
    df |> sens({{truth}}, {{estimate}}),
    df |> ppv({{truth}}, {{estimate}}),
    df |> precision({{truth}}, {{estimate}}),
    df |> recall({{truth}}, {{estimate}}),
  )
}

get_measures_simple <- function(df, truth, estimate){
  bind_rows(
    df |> accuracy({{truth}}, {{estimate}}),
    df |> spec({{truth}}, {{estimate}}),
    df |> sens({{truth}}, {{estimate}}),
    df |> precision({{truth}}, {{estimate}}),
    df |> recall({{truth}}, {{estimate}}),
  )
}


measures_by_group <- test_preds |> 
  group_by(subfamily) |> 
  get_measures_simple(truth = subfamily, estimate = predicted) |> 
  filter(!.metric == 'kap') |> 
  select(-.estimator) |> 
  filter(!is.na(.estimate)) 

measures_by_group |> 
  pivot_wider(names_from = .metric, values_from = .estimate) |> 
  View()

measures_by_group |> 
  ggplot(aes(y = fct_rev(subfamily), x = .estimate)) + 
  geom_point() + 
  facet_wrap(~.metric) + 
  labs(x = 'Score', y = 'Class',
       subtitle = 'Max-score / gathering-threshold rule-based classification measures')
## TODO decision tree: simple model
