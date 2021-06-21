library(tidyverse)
library(Biostrings)

# read in smart db & alignments for each subfamily
smart <- read_rds('./data/smart_df.rds') |> 
  select(-contains('fasta')) |> 
  mutate(align_path = Sys.glob('./data/SMART/domain_alignments/*'),
         align_fasta = map(align_path, readAAStringSet),
         align_seq = map(align_fasta, paste))
# glimpse(smart)

# unnest domains; create parent id from fasta headers
domains <- smart |> 
  select(subfamily, dom_name, dom_seq, align_seq) |>
  unnest(cols = c(dom_name, dom_seq, align_seq)) |> 
  mutate(
    dom_subfam = subfamily,
    dom_name = str_replace(dom_name, 'BPP-1', 'BPP_1'),
    dom_name = str_remove_all(dom_name, 'smart.*?-|uniprot\\|'),
    parent_id = str_remove(dom_name, '\\/.*+'),
    
    dom_id = row_number()
    ) |> 
  relocate(parent_id)

domains |> View()

length(unique(domains$dom_name))
length(unique(domains$parent_id))

# unnest proteins; create parent id from fasta headers
proteins <- smart |> 
  select(subfamily, prot_name, prot_seq) |> 
  unnest(cols = c(prot_name, prot_seq)) |> 
  mutate(
    prot_subfam = subfamily,
    prot_name = str_replace(prot_name, 'BPP-1', 'BPP_1'),
    parent_id = str_remove(prot_name, '[ ].*+'),
    prot_id = as.character(row_number())) |>
  relocate(parent_id)
proteins

# 120,549 x 4
proteins |> select(-prot_id) |>  distinct()
# 117,692 x 3
proteins |> select(-subfamily, -prot_id) |>  distinct()
 
proteins |> 
  select(-subfamily, prot_seq) |> 
  left_join(domains |> select(-c(subfamily, contains('seq')), 
                              by = c('parent_id'))) |> 
  select(-contains('seq')) |> 
  group_by(parent_id, prot_name) |> 
  mutate(prot_subs = list(prot_subfamily),
          dom_name = list(dom_name)) |> 
  mutate(len = map_int(subs, length)) |> 
  filter(len > 1) |> 
  unnest(cols = c(subs, dom_name)) |> 
  ungroup() |> 
  mutate(across(c(parent_id, prot_name), ~str_remove(.x, 'uniprot\\|'))) |> 
  relocate(subs, contains('id')) |> 
  View()
  
nrow(proteins)
length(unique(proteins$prot_name))
length(unique(proteins$prot_id))
length(unique(proteins$prot_seq))

skimr::skim(proteins)

proteins |> 
  count(parent_id) |> 
  arrange(desc(n))

proteins |> 
  dplyr::
# glimpse(domains)
# glimpse(proteins)

smart_seqs <- 
  left_join(proteins, domains, by = c('subfamily', 'parent_id')) |> 
  group_by(parent_id) |> 
  summarize(across(c(prot_name, prot_seq, dom_name, dom_seq, align_seq),
                   ~list(.x))) |> 
  ungroup()

glimpse(smart_seqs)
rm(domains, proteins, smart)




# testing <- 
#   smart_seqs |> 
#   group_by(subfamily) |> 
#   slice_sample(prop = 0.25) |> 
#   ungroup()
# 
# training <- 
#   anti_join(smart_seqs, testing, by = c("parent_id", "subfamily", "data"))
# 
# glimpse(testing)
# glimpse(training)
# 
# full_join(
#   training |> count(subfamily, name = 'training'),
#   testing |> count(subfamily, name = 'testing')
# )

# # write out an aligned fasta for a given subfamily
# write_fasta <- function(path, seqs, names){
#   seqs <- AAStringSet(seqs)
#   names(seqs) <- names
#   writeXStringSet(seqs, path)
# }
# 
# # store in /data/
# system('mkdir ./data/training_aligns')
# 
# # unnest and then re-nest individual lists
# training |> 
#   unnest(data) |> 
#   group_by(subfamily) |> 
#   mutate(path = paste0('./data/training_aligns/', subfamily, '.aln')) |> 
#   pull(path)

