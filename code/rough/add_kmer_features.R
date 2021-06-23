library(tidyverse)
library(Biostrings)
library(kmer)


# need to convert sequence from character to AAbin object
as_AAbin <- function(z, simplify = FALSE){
  res <- if(length(z) == 1 & simplify) charToRaw(z) else lapply(z, charToRaw)
  attr(res, "rerep.names") <- attr(z, "rerep.names")
  attr(res, "rerep.pointers") <- attr(z, "rerep.pointers")
  class(res) <- "AAbin"
  return(res)
}

# converts seq to AAbin obj, computes kmer proportions
add_kmers <- function(df, seq){
  kmers <- df |> 
    pull({{seq}}) |> 
    as_AAbin() |> 
    # count kmers and take proportion
    kcount(k = 2) |> 
    apply(1, function(x) x/sum(x)) |> 
    t() |> 
    as_tibble()
  df |> bind_cols(kmers)
}

# read datasets with sequences and hmmscores
integrases <- read_rds('./data/smart_refseqs_hmm_scores.rds') |> 
  select(-id) |> 
  relocate(subfamily, contains('prot')) 

non_integrases <- read_rds('./data/non_integrases_hmm_scores.rds') |> 
  mutate(prot_seq = seq,
         subfamily = factor('non_integrase'),
         prot_name = trimws(prot_name)) |> 
  select(-c(name, seq, title)) |> 
  relocate(subfamily, contains('prot'))
  
glimpse(integrases)
glimpse(non_integrases)

# join integrase and non-integrase data
# make all sequences uppercase
# remove any obs with non std AA residues in sequence
# add kmer proportions as columns
full_data <- 
  bind_rows(integrases, non_integrases) |>
  mutate(prot_seq = toupper(prot_seq)) |> 
  filter(!str_detect(prot_seq, 'B|J|O|U|X|Z|-')) |> 
  add_kmers(prot_seq)

write_rds(full_data, './data/full_classifier_dataset.rds', compress = 'gz')


