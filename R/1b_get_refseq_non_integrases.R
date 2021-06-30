### Script to retrieve non-integrases from RefSeq
### Get non-integrase proteins for 'other' class for classifier to serve as negative examples.

## Libs ---------------------------------------------------------------

library(furrr)
library(Biostrings)
library(tidyverse)
library(rentrez)
library(beepr)

## Functions for NCBI downloads ----------------------------------------------

# `get_ncbi_ids` does a basic NCBI search, retrieves metadata and record ids. Sets up the fetching process in get_ESummary_df and get_Efasta functions.
get_ncbi_ids <- function(searchexp, db){
  
  # find out how many ids available from 1st search
  Esearch <- rentrez::entrez_search(
    db = db, 
    term = searchexp, 
    retmax = 100)
  # Use the 'count' from first search to get all ids;
  # Get the web history token with 'use_history = TRUE'
  Esearch2 <- rentrez::entrez_search(
    db = db, 
    term = searchexp, 
    retmax = Esearch$count,
    use_history = TRUE)
  return(Esearch2)
}

# get_Efasta gets all sequences for a search expression,
# returns a dataframe with title and sequence columns
get_Efasta <- function(searchexp, db, apikey, fetch_max = 1600,
                       verbose = TRUE) {
  basic_search <- get_ncbi_ids(searchexp, db = db)
  webhist <- basic_search$web_history
  id.count <- min(basic_search$count, fetch_max)
  message('Search expression: ', searchexp)
  message(basic_search$count, ' records are available. Fetching ', id.count, ' fasta records')
  
  # Get 500 fasta records per iteration, concatenate them
  fasta <- furrr::future_map(
    seq(1, id.count, 500),
    ~{
      entrez_fetch(
        db = db, 
        web_history = webhist,
        retstart = .x-1, 
        retmax=500, 
        api_key = apikey,
        rettype = 'fasta')
    }) |> 
    purrr::reduce(paste0)
  
  # split fasta into lines and write a temp file
  tmp <-  tempfile()
  write(fasta, tmp)
  # read lines as a DNAStringSet object using Biostrings package
  fasta <- Biostrings::readBStringSet(tmp)
  # arrange fasta into df
  fasta.df <- tibble(title = names(fasta), seq = paste(fasta))
  
  if (verbose) 
    print(fasta)
  message('Finished query')
  
  return(list(fasta = fasta, df = fasta.df))
}



## list of proteins --------------------------------------------------------
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


## Downloads ---------------------------------------------------------------

# send queries to entrez; max n seqs returned is set to 1k
non_integrases <- 
  map(queries, 
      ~get_Efasta(.x, db = 'protein', verbose = F,
                  apikey = '889fdb9786a14019a0a1257196a09ba4ba08',
                  fetch_max = 1000)) |>
  enframe(name = "name", value = "value") |> 
  unnest_wider(value)

## Save --------------------------------------------------------------------
write_rds(non_integrases, './data/non_integrase_seqs/refseq_non_integrases_raw.rds')
