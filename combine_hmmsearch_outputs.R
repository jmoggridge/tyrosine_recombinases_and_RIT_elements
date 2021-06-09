library(tidyverse)
library(Biostrings)
library(rhmmer)

# All the full protein sequences for 20 subfamilies from SMART
# Unnest and split protein names to match hmmsearch results
integrases <- 
  read_rds('./data/SMART_integrases.rds') |> 
  select(subfamily, id, prot_name) |> 
  unnest(cols = prot_name) |> 
  mutate(temp = trimws(str_extract(prot_name, '.*? ')),
         prot_description = str_remove(prot_name,  '.*?\\s'),
         prot_name = temp) |> 
  select(-temp)

# want to join the 20 scores for each subfamily to each ref seq
integrases


# read hmmsearch tblout result
# keep only the best_domain_score column, renamed to the subfamily name
read_hmmersearch_tbl <- function(path){
  read_table(
    file = path,
    na = '-',
    comment = '#',
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

hmmsearches <- 
  # read all the results and join into one table
  Sys.glob("./data/SMART/hmmsearch_res/*.search.tbl") |> 
  map(read_hmmersearch_tbl) |> 
  purrr::reduce(left_join, by = 'prot_name') |> 
  # join to ref_integrases table
  right_join(integrases, by = 'prot_name') |> 
  select(prot_name, prot_description, subfamily, id, everything())

glimpse(hmmsearches)


hmmsearches |> 
  pivot_longer(cols = c(Arch1:Xer), 
               names_to = 'hmm_dom', values_to = 'hmm_score') |> 
  mutate(hmm_score = replace_na(hmm_score, 0)) |> 
  mutate(highlight = ifelse(hmm_dom == subfamily, T, F)) |> 
  ggplot(aes(x = hmm_score, y = fct_rev(hmm_dom),
             color = highlight)) +
  geom_violin() +
  facet_wrap(~subfamily) +
  labs(x = 'Best domain bit score', y = 'Hidden Markov model', 
       subtitle = 'HMM bit score profiles for each set of tyrosine integrases from the SMART subfamilies') +
  theme_light() + 
  theme(legend.position = 'null', 
        axis.text.y = element_text(size = 5))


# # all hmmsearch output
# hmmsearch_rs <- 
#   tibble(path = Sys.glob('./data/SMART/hmmsearch_res/*.tbl')) |> 
#   mutate(domain_name = str_extract(path, '[^\\s|!/]+.search.tbl'),
#          domain_name = str_remove(domain_name, '.search.tbl')) |> 
#   select(domain_name, path) |> 
#   arrange(domain_name) |> 
#   mutate(result = map(path, ~read_hmmertbl(.x)))
#   # unnest(cols = c(result))
# 
# glimpse(hmmsearch_rs)
# 
# 
# search_tbl <-
#   rhmmer::read_tblout("./data/SMART/hmmsearch_res/Arch1.search.tbl") |> 
#   select(-contains('accession'), -contains('domain_number_')) |>
#   mutate(protein = paste(domain_name, description))
# 
# names(search_tbl)
# glimpse(search_tbl)
# 

