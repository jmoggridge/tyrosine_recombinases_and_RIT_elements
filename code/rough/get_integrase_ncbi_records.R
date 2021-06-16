library(tidyverse)
library(rentrez)

apikey <- "889fdb9786a14019a0a1257196a09ba4ba08"

get_ncbi_ids <- function(searchexp, db){
  
  # find out how many ids available for search
  Esearch <- entrez_search(
    db = db, 
    term = searchexp, 
    retmax = 100
    )
  # uses 'count' value from search to get all ids;
  # get webenv token for downloading the records with 'use_history = TRUE'
  Esearch2 <- entrez_search(
    db = db, 
    term = searchexp, 
    retmax = Esearch$count,
    use_history = TRUE
    )
  return(Esearch2)
}

search_expr <- "(tyrosine[All Fields] AND recombinase[All Fields]) AND alive[prop]"
basic_search <- get_ncbi_ids(search_expr, db = 'gene')

# basic_search$count
# basic_search$ids


get_ESummary_df <- function(searchexp, db, apikey){
  basic_search <- get_ncbi_ids(searchexp, db=db)
  webhist <- basic_search$web_history
  id.count <- basic_search$count
  df <- tibble(id = basic_search$ids)
  # message('Getting summaries....')
  Esummary_list <- get_Esummaries(webhist, id.count, apikey)
  # message('done downloads....\nflattening lists....')
  df <- df %>%
    mutate(listcol = Esummary_list) %>%
    unnest_wider(listcol)
  return(df)
}

# df <- get_ESummary_df(searchexp = )


