
## functions I wrote get ncbi fasta records
source('./code/ncbi_entrez_functions.R')

# set api key and choose database
apikey <- '889fdb9786a14019a0a1257196a09ba4ba08'
db <- 'protein'

queries <- list(
  'Serine recombinase' = "serine recombinase",
  'Holliday junction resolvase' = 'Holliday junction resolvase',
  'Topoisomerase' = 'topoisomerase',
  'DDE transposase' = "DDE transposase",
  'DNA helicase recQ' = 'DNA helicase recQ',
  'Thermonuclease' = 'thermonuclease',
  'Ribonuclease' = 'ribonuclease',
  'Helicase' = 'helicase',
  'Histone' = 'histone',
  'Argonaute' = 'argonaute',
  'Exodeoxyribonuclease' = 'Exodeoxyribonuclease',
  'Restriction enz' = 'restriction endonuclease',
  'MbeA' = 'DNA relaxase mbeA domain protein',
  'RpnA' = 'Recombination-promoting nuclease RpnA',
  'RpnB' = 'Recombination-promoting nuclease RpnB',
  'RpnC' = 'recombination-promoting nuclease RpnC',
  'GamL' = 'host nuclease inhibitor GamL',
  'Nuclease SbcCD subunit D' = 'Nuclease SbcCD subunit D',
  'NikB' = 'NikB',
  'RecB' = 'exodeoxyribonuclease V subunit beta',
  'RecC' = 'exodeoxyribonuclease V subunit gamma',
  'RecD' = 'exodeoxyribonuclease V subunit alpha',
  'RecJ' = 'single-stranded-DNA-specific exonuclease RecJ',
  'RecN' = 'DNA repair protein RecN',
  'RecR' = 'recombination protein RecR',
  'Exonuclease V' = 'Exonuclease V',
  'IS607' = 'IS607 family transposase',
  'DnaX' = 'DNA polymerase III subunit gamma/tau',
  'Gntr' = 'gntr family transcriptional regulator',
  'Adenylosuccinate lyase' = 'Adenylosuccinate lyase',
  'IhfA' = 'integration host factor subunit alpha',
  'IhfB' = 'Integration host factor subunit beta',
  'radA' = 'DNA repair protein RadA'
  
) |> 
  map(~paste0('(', .x, '[Protein Name]) AND refseq[filter]'))
queries


other_seqs <- 
  map(queries, ~ get_Efasta(.x, db, apikey, fetch_max = 1000)) |> 
  enframe(name = "name", value = "value") |> 
  unnest_wider(value)

glimpse(other_seqs)

other_seqs |> 
  unnest(df)

write_rds(other_seqs, './data/non_integrases_df.rds')


# 
# queries2 <- list(
#   'RecR' = 'recombination protein RecR',
#   'RecN' = 'DNA repair protein RecN'
#   ) |> 
#   map(~paste0('(', .x, '[Protein Name]) AND refseq[filter]'))
# 
# other_seqs2 <- 
#   map(queries2, ~ get_Efasta(.x, db, apikey, fetch_max = 1000)) |> 
#   enframe(name = "name", value = "value") |> 
#   unnest_wider(value)





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

