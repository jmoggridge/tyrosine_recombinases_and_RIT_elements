library(tidyverse)
library(Biostrings)
library(kmer)

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


# join integrase and non-integrase
# remove any obs with non std AA residues in sequence
full_data <- 
  bind_rows(integrases, non_integrases) |>
  mutate(prot_seq = toupper(prot_seq)) |> 
  filter(!str_detect(prot_seq, 'B|J|O|U|X|Z|-')) 

# smallest class size
smallest <- full_data |> group_by(subfamily) |> count() |> pull(n) |> min()

# downsample to size of smallest class
set.seed(1)

df <- full_data |> 
  group_by(subfamily) |> 
  sample_n(size =  smallest, replace = F) |> 
  ungroup()

glimpse(df)
rm(smallest, integrases, non_integrases, full_data)



# add kmer profiles

# need to convert sequence from character to AAbin object
as_AAbin <- function(z, simplify = FALSE){
  res <- if(length(z) == 1 & simplify) charToRaw(z) else lapply(z, charToRaw)
  attr(res, "rerep.names") <- attr(z, "rerep.names")
  attr(res, "rerep.pointers") <- attr(z, "rerep.pointers")
  class(res) <- "AAbin"
  return(res)
}

seqs <- df$prot_seq
names(seqs) <- df$prot_name
aabin <- seqs |> as_AAbin()
kc <- kcount(aabin, k = 2) 
dim(kc)
kc2 <- t(apply(kc, 1, function(x) x/sum(x)))
dim(kc2)
kc3 <- kc2 |> 
  as.data.frame() |> 
  as_tibble(prot_name = rownames())
kc3

kmer_count <- 
  df$prot_seq |> 
  as_AAbin() |> 
  # count kmers and take proportion
  kcount(k = 2) |> 
  apply(1, function(x) x/sum(x)) |> 
  t() |> 
  as_tibble()

df <- df |> 
  bind_cols(kmer_count)
df

