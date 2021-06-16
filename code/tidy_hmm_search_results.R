## Tidy hmmsearch results for integrases and non-integrases, then prepare data for classifier.
# I parsed all the `hmmsearch` results tables and rectangled the data by joining the bitscores to the protein names with their true subfamilies in the `SMART_db` dataframe.  
# In the second part, I do the same with the non-integrase examples

library(tidyverse)
library(purrr)

## fx ----
read_hmmersearch_tbl <- function(path){
  # read hmmsearch tblout result
  # keep only the best_domain_score column, 
  # rename it as the subfamily name for the hmm
  read_table(
    file = path,
    na = '-',
    comment = '#',
    col_types = cols(),
    col_names =  c(
      "prot_name", "prot_accession", 
      "hmm_name", "hmm_accesion",
      "seq_eval", "seq_score", "seq_bias",
      "best_dom_eval", "best_dom_score", "best_dom_bias",
      'domain_n_exp', 'domain_n_reg', 'domain_n_clu', 'domain_n_ov',
      'domain_n_env', 'domain_n_dom', 'domain_n_rep', 'domain_n_inc',
      "description")
  ) |> 
    select(prot_name, hmm_name, best_dom_score) |> 
    distinct() |> 
    pivot_wider(names_from = hmm_name, values_from = best_dom_score) |> 
    unnest(cols = c())
}

## integrases ----

# All the full protein sequences for 20 subfamilies from SMART
# Unnest and split protein names to match hmmsearch results
integrases <- 
  read_rds('./data/SMART_db.rds') |> 
  select(subfamily, id, prot_name, prot_seq) |> 
  unnest(cols = c(prot_name, prot_seq)) |> 
  mutate(temp = trimws(str_extract(prot_name, '.*? ')),
         prot_description = str_remove(prot_name,  '.*?\\s'),
         prot_name = temp) |> 
  select(-temp)

glimpse(integrases)

# read all the results and join into one table
# join to ref_integrases table by protein names, such that each row has the
# scores across all models for a given sequence
integrases <- 
  Sys.glob("./data/SMART/hmmsearch_res/*.search.tbl") |> 
  map(read_hmmersearch_tbl) |> 
  purrr::reduce(full_join, by = 'prot_name') |> 
  right_join(integrases, by = 'prot_name') |> 
  select(prot_name, prot_description, subfamily, id, everything())

glimpse(integrases)

write_rds(integrases, './data/smart_refseqs_hmm_scores.rds', compress = 'gz')


## plots -----
# I plotted the score profiles for each subfamily to explore the discriminative power of the bit-scores. This shows the distribution of scores for a set of sequences from a given subfamily (panel) against each subfamily hmm (y-axis), in terms of the bit scores (x-axis). Mostly the true subfamily has the best-scoring HMM alignment but in the case of `Int_Tn916` and `Int_BPP-1` integrases, the scores are very similar between the two HMMs.

plt1 <- integrases |> 
  pivot_longer(cols = c(Arch1:Xer), 
               names_to = 'hmm_dom', values_to = 'hmm_score') |> 
  mutate(hmm_score = replace_na(hmm_score, 0)) |> 
  mutate(highlight = ifelse(hmm_dom == subfamily, T, F)) |> 
  ggplot(aes(x = hmm_score, y = fct_rev(hmm_dom),
             color = highlight)) +
  geom_violin() +
  labs(x = 'Best domain bit score', y = 'Hidden Markov model', 
       subtitle = 'HMM bit score profiles for each set of tyrosine integrases from the SMART subfamilies') +
  facet_wrap(~subfamily, nrow = 3) +
  theme_light() + 
  theme(legend.position = 'null', 
        axis.text.y = element_text(size = 5),
        panel.grid.major.x = element_blank()
  )

# We could also exchange the HMMs and true subfamilies in the plot above to get the scores for each model as a panel, with the true classes on the y-axis.

plt2 <- integrases |> 
  pivot_longer(cols = c(Arch1:Xer), 
               names_to = 'hmm_dom', values_to = 'hmm_score') |> 
  mutate(hmm_score = replace_na(hmm_score, 0)) |> 
  mutate(highlight = ifelse(hmm_dom == subfamily, T, F)) |> 
  ggplot(aes(x = hmm_score, y = fct_rev(subfamily),
             color = highlight)) +
  geom_violin() +
  # ggthemes::geom_tufteboxplot(stat = "fivenumber",
  #                             size = 1.5, whisker.type = 'line') +
  facet_wrap(~hmm_dom, nrow = 2) +
  labs(x = 'Best domain bit score', y = 'True subfamily', 
       subtitle = 'HMM bit score profiles for each hidden Markov model') +
  theme_light() + 
  theme(legend.position = 'null', 
        axis.text.y = element_text(size = 5))





## Non-integrase seqs ----

# had to alter the parsing a bit for these
read_hmmersearch_tbl2 <- function(path){
  read_delim(
    file = path,
    na = '-',
    delim = '  ',
    trim_ws = T,
    comment = '#',
    col_types = cols(),
    col_names =  c(
      "prot_name", "prot_accession", 
      "hmm_name", "hmm_accesion",
      "seq_eval", "seq_score", "seq_bias",
      "best_dom_eval", "best_dom_score", "best_dom_bias",
      'domain_n_exp', 'domain_n_reg', 'domain_n_clu', 'domain_n_ov',
      'domain_n_env', 'domain_n_dom', 'domain_n_rep', 'domain_n_inc',
      "description")
  ) |> 
    select(prot_name, hmm_name, best_dom_score) |> 
    distinct() |> 
    pivot_wider(names_from = hmm_name, values_from = best_dom_score) |> 
    unnest(cols = c())
}

# ./data/non_integrase_seqs/hmmsearch_res/
non_int_rs <- 
  Sys.glob('./data/non_integrase_seqs/hmmsearch_res/*.tbl') |> 
  map(read_hmmersearch_tbl2) |> 
  purrr::reduce(full_join, by = 'prot_name')

non_integrases <- 
  # read non-int data
  read_rds('./data/non_integrases_df.rds') |> 
  # WP_161670389.1 has an integrase domain, discard
  filter(!prot_name == 'WP_161670389.1') |> 
  left_join(y = non_int_rs, 
            by = 'prot_name') |> 
  mutate(across(Arch1:Xer, ~replace_na(.x, 0)))

# hmm scores of non_integrases are very low indeed.
non_integrases |> 
  summarise(across(.cols = Arch1:Xer, ~mean(.x))) |> 
  pivot_longer(everything())
non_integrases |> 
  summarise(across(.cols = Arch1:Xer, ~max(.x))) |> 
  pivot_longer(everything())

non_integrases |> 
  filter(`Int_BPP-1` == max(`Int_BPP-1`)) |> 
  pull(prot_name)
non_integrases |> 
  filter(`Int_Tn916` == max(`Int_Tn916`)) |>
  pull(prot_name)


write_rds(non_integrases, './data/non_integrases_hmm_scores.rds', compress = 'gz')
