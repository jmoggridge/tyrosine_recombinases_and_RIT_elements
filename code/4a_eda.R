library(tidyverse)
library(furrr)

# TODO EDA of prepped data

test <- read_rds('./data/test_df_scored.rds')
train <- read_rds('./data/train_df_scored.rds')
both <- bind_rows(train |> mutate(data = 'train'),
                  test |> mutate(data  = 'test')) |> 
  select(subfamily, prot_seq, data)



# check for any weird characters: XBZJ
map(list(train, test), 
    ~{.x |> 
        select(subfamily, acc, prot_seq) |> 
        mutate(n_nonAA = str_count(prot_seq, '[^ACDEFGHIKLMNPQRSTVWYX]'),
               nonAA = str_extract_all(prot_seq, '[^ACDEFGHIKLMNPQRSTVWYX]')) |> 
        unnest(cols = c(nonAA)) |> 
        count(nonAA)
    }
)


## plot seq lengths in test + train data, by subfamily (or non_integrase)
both |> 
  mutate(seqlen = nchar(prot_seq)) |> 
  ggplot(aes(x = seqlen, y = fct_rev(subfamily))) +
  geom_jitter(alpha = 0.05, size = 0.05) +
  scale_x_log10() +
  facet_wrap(~data) +
  labs(x = 'residues', y ='', subtitle = 'Protein sequence length by class') +
  theme_classic()
  
aa <- list('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','X')

# get AA profiles
count_aa <- function(seq){
  tibble(aa = aa) |> 
    mutate(count = map_dbl(aa, ~str_count(seq, .x)) / nchar(seq)) |> 
    pivot_wider(names_from = aa, values_from = count)
}

plan(multisession, workers = availableCores())
aa_pro <- both |>
  mutate(profile = future_map(prot_seq, count_aa)) |>
  unnest(profile) |> 
  mutate_if(is.numeric, ~ {(.x - mean(.x)) / sd(.x)} )

# plot AA profiles
aa_pro |> 
  filter(!subfamily == 'non_integrase') |> 
  pivot_longer(A:Y, names_to = 'AA', values_to = 'proportion') |> 
  ggplot(aes(AA, proportion)) +
  geom_violin(scale = 'width', draw_quantiles = c(0.5), trim = T) +
  facet_wrap(~subfamily) +
  labs(x = 'residue', y = 'normalized proportion (vs all seqs)')

# plot AA profiles
aa_pro |> 
  mutate(group = ifelse(subfamily == 'non_integrase', 'non_int', 'int')) |> 
  pivot_longer(A:Y, names_to = 'AA', values_to = 'proportion') |> 
  ggplot(aes(AA, proportion)) +
  geom_violin(scale = 'width', draw_quantiles = c(0.5), trim = T) +
  facet_wrap(~group) +
  labs(x = 'residue', y = 'normalized proportion (vs all seqs)')


# hmmsearch score profiles: each panel is the results for a given subfamily's sequences
bind_rows(
  train |> mutate(data = 'train'),
  test |> mutate(data  = 'test')
  ) |> 
  select(subfamily, Arch1:Xer) |> 
  pivot_longer(cols = c(Arch1:Xer), 
               names_to = 'hmm_dom', values_to = 'hmm_score') |> 
  mutate(hmm_score = replace_na(hmm_score, 0)) |> 
  mutate(highlight = ifelse(hmm_dom == subfamily, T, F)) |> 
  ggplot(aes(x = hmm_score, y = fct_rev(hmm_dom),
             color = highlight)) +
  # geom_violin() +
  geom_boxplot(outlier.size = 0.2) +
  facet_wrap(~subfamily, nrow = 3) +
  labs(x = 'Best domain bit score', y = 'Hidden Markov model', 
       subtitle = 'HMM bit score profiles for each set of tyrosine integrases from the SMART subfamilies') +
  theme_light() + 
  theme(legend.position = 'null', 
        axis.text.y = element_text(size = 5),
        panel.grid.major.x = element_blank()
        )




# plot gathering thresholds on pointrange of self-scores
raw_self_scores <- train |> 
  select(subfamily, Arch1:Xer) |> 
  pivot_longer(cols = Arch1:Xer, names_to = 'hmm_name', values_to = 'hmm_score') |> 
  filter(subfamily == hmm_name) |> 
  ggplot(aes(x = fct_rev(hmm_name), y = hmm_score)) +
  geom_boxplot(outlier.size = 0.2) +
  stat_summary(fun = min, geom = 'crossbar', color ='red', size = 0.35) +
  labs(x ='model', y = 'self-scores (raw bit scores)', 
       subtitle = 'HMM training sequences self-scoring and gathering thresholds') +
  theme_classic() +
  coord_flip() 


# normalized gathering thresholds and 
normalize_self_scores <- train |> 
  select(subfamily, Arch1:Xer) |>
  mutate(across(Arch1:Xer, ~c(scale(.x, center = T, scale = T)), .names = '{col}')) |>
  pivot_longer(cols = Arch1:Xer, names_to = 'hmm_name', values_to = 'hmm_score') |> 
  filter(subfamily == hmm_name) |> 
  
  ggplot(aes(x = fct_rev(hmm_name), y = hmm_score)) +
  geom_boxplot(outlier.size = 0.2) +
  stat_summary(fun = min, geom = 'crossbar', color ='red', size = 0.35) +
  labs(x ='model', y = 'self-scores (normalized vs all scores for each model)', 
       subtitle = 'HMM training sequences self-scoring and gathering thresholds') +
  theme_classic() +
  coord_flip() 

library(patchwork)
raw_self_scores + normalize_self_scores