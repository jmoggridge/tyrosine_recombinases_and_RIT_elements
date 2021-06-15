
## get ncbi fasta records

source('./code/ncbi_entrez_functions.R')

# set api key and choose database
apikey <- '889fdb9786a14019a0a1257196a09ba4ba08'
db <- 'protein'

queries <- list(
  'thermonuclease' = 'thermonuclease',
  'ribonuclease' = 'ribonuclease',
  'histone' = 'histone',
  'argonaute' = 'argonaute',
  'restriction enz' = 'restriction endonuclease',
  'MbeA' = 'DNA relaxase mbeA domain protein',
  'Exodeoxyribonuclease' = 'Exodeoxyribonuclease',
  'RpnB' = 'Recombination-promoting nuclease RpnB',
  'RpnC' = 'recombination-promoting nuclease RpnC',
  'GamL' = 'host nuclease inhibitor GamL',
  'Nuclease SbcCD subunit D' = 'Nuclease SbcCD subunit D',
  'helicase' = 'helicase',
  'NikB' = 'NikB',
  'RecB' = 'RecB',
  'RecC' = 'RecC',
  'RecD' = 'RecD',
  'Holliday junction resolvase' = 'Holliday junction resolvase',
  'Gntr' = 'gntr family transcriptional regulator',
  'topoisomerase' = 'topoisomerase',
  'DDE transposase' = "DDE transposase",
  'DNA helicase recQ' = 'DNA helicase recQ',
  'Adenylosuccinate lyase' = 'Adenylosuccinate lyase',
  'serine recombinase' = "serine recombinase"
) |> 
  map(~paste0('(', .x, '[Protein Name]) AND refseq[filter]'))
queries

queries2 <- list(
  'RecC' = 'exodeoxyribonuclease V subunit gamma',
  'RecD' = 'exodeoxyribonuclease V subunit alpha',
  'RecB' = 'exodeoxyribonuclease V subunit beta',
  'Exonuclease V' = 'Exonuclease V',
  'helicase AddB' = 'helicase AddB',
  'dnaX' = 'DNA polymerase III subunit gamma/tau'
  
  ) |> 
  map(~paste0('(', .x, '[Protein Name]) AND refseq[filter]'))

other_seqs2 <- 
  map(queries2, ~ get_Efasta(.x, db, apikey, fetch_max = 1000)) |> 
  enframe(name = "name", value = "value") |> 
  unnest_wider(value)


other_seqs <- 
  map(queries, ~ get_Efasta(.x, db, apikey, fetch_max = 1000)) |> 
  enframe(name = "name", value = "value") |> 
  unnest_wider(value)

glimpse(other_seqs)




# downloaded all the relaxase sequences from pfam

seqs <- readAAStringSet('./data/non_integrase_sequences/relaxase_PF03432_full.fa.txt')
seqs_df <- tibble(
  name = 'relaxase',
  fasta = list(seqs)) |> 
  nest(fasta = fasta) |> 
  mutate(
  df = list(tibble(title = names(seqs),
              seq = paste0(seqs))))

glimpse(seqs_df)
seqs |> unnest_wider(col = df)
# should remove this sequence: A0A1K1NN93_9FIRM, has a phage integrase domain

