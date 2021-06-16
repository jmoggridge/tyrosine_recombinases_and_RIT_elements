### Script to retrieve non-integrases from RefSeq

library(furrr)
library(Biostrings)
library(tidyverse)

# functions I wrote get ncbi fasta records
source('./code/ncbi_entrez_functions.R')

# list of proteins
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
  'Restriction enz' = 'restriction endonuclease',
  'MbeA' = 'DNA relaxase mbeA domain protein',
  'RpnA' = 'Recombination-promoting nuclease RpnA',
  'RpnB' = 'Recombination-promoting nuclease RpnB',
  'RpnC' = 'recombination-promoting nuclease RpnC',
  'GamL' = 'host nuclease inhibitor GamL',
  'Nuclease SbcCD subunit D' = 'Nuclease SbcCD subunit D',
  'NikB' = 'NikB',
  'Exodeoxyribonuclease' = 'Exodeoxyribonuclease',
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
  'RadA' = 'DNA repair protein RadA'
) |> 
  # add search syntax
  map(~paste0('(', .x, '[Protein Name]) AND refseq[filter]'))

# send queries to entrez; max n seqs returned is set to 1k
non_integrases <- 
  map(queries, 
      ~get_Efasta(.x, db = 'protein', verbose = F,
                  apikey = '889fdb9786a14019a0a1257196a09ba4ba08',
                  fetch_max = 1000)) |>
  enframe(name = "name", value = "value") |> 
  unnest_wider(value)

write_rds(non_integrases, './data/non_integrase_seqs/refseq_non_integrases_raw.rds')
