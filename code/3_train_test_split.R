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
    dom_name = str_replace(dom_name, 'BPP-1', 'BPP_1'),
    parent_id = str_remove(dom_name, 'smart.*?-') |> str_remove('\\/.*+')
    ) |> 
  relocate(parent_id)

# unnest proteins; create parent id from fasta headers
proteins <- smart |> 
  select(subfamily, prot_name, prot_seq) |> 
  unnest(cols = c(prot_name, prot_seq)) |> 
  mutate(prot_name = str_replace(prot_name, 'BPP-1', 'BPP_1'),
         parent_id = str_remove(prot_name, '[ ].*+')) |> 
  relocate(parent_id)

# glimpse(domains)
# glimpse(proteins)

smart_seqs <- 
  left_join(proteins, domains, by = c('subfamily', 'parent_id')) |> 
  group_by(parent_id) |> 
  nest(data = c(prot_name, prot_seq, dom_name, dom_seq, align_seq)) |> 
  ungroup()

glimpse(smart_seqs)
rm(domains, proteins, smart)

testing <- 
  smart_seqs |> 
  group_by(subfamily) |> 
  slice_sample(prop = 0.25) |> 
  ungroup()

training <- 
  anti_join(smart_seqs, testing, by = c("parent_id", "subfamily", "data"))

glimpse(testing)
glimpse(training)

full_join(
  training |> count(subfamily, name = 'training'),
  testing |> count(subfamily, name = 'testing')
)

# write out aligned fasta for hmms
write_fasta <- function(path, seqs, names){
  seqs <- AAStringSet(seqs)
  names(seqs) <- names
  writeXStringSet(seqs, path)
}



