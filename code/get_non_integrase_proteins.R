### Get non-integrase proteins for 'other' class for classifier

library(tidyverse)
library(furrr)
library(rentrez)
library(beepr)
library(Biostrings)

#### Functions for NCBI downloads ----

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
  message('Returning metadata from basic search...')
  print(Esearch2)
  return(Esearch2)
}
# get_ESummary_df downloads all NCBI nucleotide summary records for a 
# search expression and returns a dataframe
get_ESummary_df <- function(searchexp, db, apikey){
  
  # perform basic search to get webhistory and records count
  basic_search <- get_ncbi_ids(searchexp, db)
  web_history <- basic_search$web_history
  id.count <- basic_search$count
  
  # init list to gather downloads 
  Esummary_list <- list()
  # init df to compile records from all downloads
  df <- tibble(id = basic_search$ids)
  
  # display downloads progress
  message('Getting summaries....')
  progress_bar = txtProgressBar(min=0, max=id.count, style = 1, char="=")
  
  # iterate until all summary records obtained
  for (i in seq(1, id.count, 500)) {
    # add each query to growing Esummary list
    Esummary_list <-  c(
      Esummary_list,
      entrez_summary(db = db, web_history = web_history,
                     retstart = i-1, retmax=500, api_key = apikey,
                     always_return_list = TRUE, retmode = 'json')
    )
    setTxtProgressBar(progress_bar, value = i)
    Sys.sleep(0.11)
  }
  message('\ndone downloads....\nflattening lists....')
  df <- df |>
    mutate(listcol = Esummary_list) |>
    unnest_wider(listcol)
  return(df)
}

# get_Efasta gets all sequences for a search expression,
# returns a dataframe with title and sequence columns
get_Efasta <- function(searchexp, db, apikey, fetch_max = 1600,
                       verbose = TRUE) {
  basic_search <- get_ncbi_ids(searchexp, db = db)
  webhist <- basic_search$web_history
  id.count <- min(basic_search$count, fetch_max)
  message(paste('Fetching', id.count, 'fasta records'))
  progress_bar = txtProgressBar(min=0, max=id.count, style = 1, char="=")
  
  # Get 500 fasta records per iteration, concatenate them
  fasta <- furrr::future_map(
    seq(1, id.count, 500),
    ~{
      setTxtProgressBar(progress_bar, value = .x)
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

#### Downloading the NCBI Apicomplexa 18S sequence dataset ----

searchexp <- "(topoisomerase[Protein Name]) AND refseq[filter]"
# set api key and choose database
apikey <- '889fdb9786a14019a0a1257196a09ba4ba08'
db <- 'protein'

queries <- list(
  'CRE recombinase' = 'CRE recombinase',
  'FLP recombinase' = 'FLP recombinase',
  'topoisomerase' = 'topoisomerase',
  'DDE transposase' = "DDE transposase",
  'serine recombinase' = "serine recombinase"
) |> 
  map(~paste0('(', .x, '[Protein Name]) AND refseq[filter]'))
queries

other_seqs <- 
  map(queries, ~ get_Efasta(.x, db, apikey, fetch_max = 1000)) 

other_seqs |> 
  enframe(name = "name", value = "value") |> 
  unnest_wider(value)

glimpse(other_seqs)
