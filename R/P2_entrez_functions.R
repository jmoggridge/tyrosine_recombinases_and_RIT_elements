library(tidyverse)
library(rentrez)
library(glue)
library(progress)
library(tictoc)
library(beepr)

set_entrez_key('889fdb9786a14019a0a1257196a09ba4ba08')
Sys.getenv('ENTREZ_KEY')

# FUNCTIONS --------------------------------------------------

## Sequence parsing ------------------------------------------

# fasta to DNA set
fasta_to_DNAss <- function(x){
  temp <- tempfile()
  write_file(x, temp)
  Biostrings::readDNAStringSet(temp)
}

# protein to AAss
fasta_to_AAss <- function(x){
  temp <- tempfile()
  write_file(x, temp)
  Biostrings::readAAStringSet(temp)
}

parse_nuc_data <- function(nuc_data){
  nuc_data |> 
    mutate(
      nuc_ss = map(get, fasta_to_DNAss),
      nuc_name = map(nuc_ss, names),
      nuc_seq = map(nuc_ss, paste)
    ) |> 
    select(id, nuc_name, nuc_seq)
}

## Link IDS across database ----------------------------------

# for a set of search terms, return cdd identifiers and summaries
# data needs to a df, term is the column of search expressions
cdd_search <- function(data, term){
  data |> 
    mutate(
      search = 
        map({{term}}, ~entrez_search(db = 'cdd', term = .x,
                                     use_history = T)),
      cdd_id = map(search, pluck, 'ids'),
    ) |> 
    unnest(cdd_id) |>
    # don't want the superfamily DNA_BRE_C id '412227'
    filter(cdd_id != '412227') |>
    mutate(
      cdd_summary = map(cdd_id, ~entrez_summary(db = 'cdd', id = .x)),
    ) |> 
    select(-search)
}

# for a set of CDD ids, returns the linked protein ids.
link_cdd_protein <- function(data, id) {
  data |>
    mutate(link = map(
      .x = {{id}},
      ~entrez_link(
        id = .x,
        dbfrom = "cdd",
        db = "protein"
      )
    ),
    prot_id = map(link, ~.x$links$cdd_protein_specific_2)
    ) |> 
    select(-link)
}

# linking protein to nucleotide data
link_protein_nuccore <- function(data, id){
  df <- data |>
    # can get db elinks from proteins to  nuccore records
    mutate(link = map(
      .x = {{id}},
      .f = ~ {
        entrez_link(
          id = .x,
          dbfrom = "protein",
          db = "nuccore"
        )
      }
    ),
    nuc_id = map(link, ~.x$links$protein_nuccore)
    ) |>
    # drop any protein records without a linked nucleotide
    mutate(drop = map_lgl(nuc_id, ~{!is.null(.x)})) |> 
    filter(drop == T) |> 
    select(-c(drop, link))
  return(df)
}

link_protein_nuccore_large <- function(data, id, chunk_size=100){
  prot_ids <- data |>
    select({{prot_id}}) |>
    distinct()
  
  n <- nrow(prot_ids)
  pb <- progress_bar$new(
    format = "  downloading [:bar] :percent eta: :eta",
    total = n/100, clear = FALSE, width= 60
  )
  
  # split prot ids into in chunks
  prot_nuc_ids <-
    tibble(
      df = map(
        .x = seq(1, n, chunk_size),
        .f = ~prot_ids |> slice(.x:min(.x + chunk_size, n))
      )
    ) |> 
    # get nuc ids
    mutate(
      nuc_id = map(df, ~{
        pb$tick()
        link_protein_nuccore(.x, id = prot_id)
      }))
  return(prot_nuc_ids)
}


# get linked taxonomy records for each nucleotide
# input tibble with a column of ids
link_nuccore_taxonomy <- function(data, id){
  df <- data |>
    select({{id}}) |>
    distinct() 
  pb <- progress_bar$new(total = nrow(df))
  df |>
    mutate(
      tax_link = map(
        .x = {{id}},
        .f = ~{
          pb$tick()
          entrez_link(
            id = .x,
            dbfrom = 'nuccore',
            db = 'taxonomy',
          )}
      ),
      tax_id = map(tax_link, ~{.x$links$nuccore_taxonomy})
    ) |>
    select(-tax_link) |>
    unnest(tax_id, keep_empty = T) |>
    right_join(data)
}



## DATA GETTERS -----------------------------------------------

## fetch  entrez data for a set of ids for a given database
fetch_data <- function(data, id, db = 'nuccore', chunk_size = 50){
  tic()
  # only keep unique identifiers
  df <- data |> 
    transmute(id = {{id}}) |> 
    distinct() 
  message('input data')
  glimpse(df)
  # split data into chunks of chunk_size
  n <- chunk_size
  dfs <- 
    tibble(
      id_set = map(.x = seq(1, nrow(df), n),
                   .f = ~df |> slice(.x:min(.x + n - 1, nrow(df))))
    ) |> 
    unnest_wider(id_set)
  
  message('Posting requests')
  pb <- progress_bar$new(total = nrow(dfs))
  dfs <- dfs |> 
    mutate(
      token = map(id, ~{
        Sys.sleep(0.2)
        pb$tick()
        entrez_post(db = db, id =.x)
      }))
  
  message('Fetching sequences')
  pb <- progress_bar$new(total = nrow(dfs))
  dfs <- dfs |> 
    mutate(
      get = map(token, ~{
        Sys.sleep(0.2)
        pb$tick()
        entrez_fetch(db = db, web_history = .x, rettype = 'fasta')})
    )
  
  message('Cleaning up')
  if (db == 'nuccore'){
    return(dfs)
  } else if (db == 'protein'){
    dfs <- dfs |>
      mutate(
        prot_ss = map(get, fasta_to_AAss),
        prot_name = map(prot_ss, names),
        prot_seq = map(prot_ss, paste)
      )}
  toc()
  return(dfs |> select(id, contains('name'), contains('seq')))
}

fetch_summaries <- function(data, id, db = 'nuccore', chunk_size = 50){
  
  tic()
  # only keep unique identifiers
  df <- data |> 
    transmute(id = {{id}}) |> 
    distinct() 
  message('input data')
  glimpse(df)
  # split data into chunks of chunk_size
  n <- chunk_size
  dfs <- 
    tibble(
      id_set = map(.x = seq(1, nrow(df), n),
                   .f = ~df |> slice(.x:min(.x + n - 1, nrow(df))))
    ) |> 
    unnest_wider(id_set)
  
  message('Posting requests')
  pb <- progress_bar$new(total = nrow(dfs))
  dfs <- dfs |> 
    mutate(
      token = map(id, ~{
        Sys.sleep(0.2)
        pb$tick()
        entrez_post(db = db, id =.x)
      }))
  
  message('Fetching summaries')
  pb <- progress_bar$new(total = nrow(dfs))
  dfs <- dfs |> 
    mutate(
      summary = map(token, ~{
        Sys.sleep(0.2)
        pb$tick()
        entrez_summary(db = db, web_history = .x)  |> 
          enframe() |> 
          unnest_wider(value)})
    )
  toc()
  return(dfs)
}

### better way to create files?
fetch_nuc_datasets <- function(data, set_number) {
  message(glue('set number: {set_number}'))
  data |> 
    fetch_data(id = nuc_id, db = 'nuccore', chunk_size = 100) |> 
    parse_nuc_data() |> 
    write_rds(nuc_data12, glue('./data/CDD/nuc_data/{set_number}.rds'))
  return('./data/CDD/nuc_data/{set_number}.rds')
}


# return taxonomic records for ids
# input dataframe with column called tax_id,
fetch_taxonomy <- function(data, id){
  pb <- progress_bar$new(total = nrow(data))
  data %>% 
    # TODO check that map_chr works here
    mutate(
      xml = 
        map_chr({{id}}, 
                ~{
                  pb$tick()
                  entrez_fetch(db = 'taxonomy', id = .x, 
                               rettype="xml")}
        ),
      tax_name = 
        str_extract(xml, '<ScientificName>.*?</ScientificName>') %>%
        str_remove_all('<ScientificName>|</ScientificName>'),
      tax_parent = 
        str_extract(xml, '<ParentTaxId>[0-9]+') %>% 
        str_remove_all('<ParentTaxId>'),
      tax_lineage = 
        str_extract(xml, '<Lineage>.*?</Lineage>\n') %>% 
        str_remove_all('<Lineage>|</Lineage>\n'),
      tax_genetic_code =
        str_extract(xml, '<GeneticCode>\n        <GCId>[0-9]+') %>% 
        str_remove_all('<GeneticCode>\n        <GCId>'),
    )  %>% 
    separate(
      tax_lineage, sep = '; ', remove = T, extra = 'merge',
      into = c(
        'no rank',
        'superkingdom',
        'clade',
        'phylum',
        'class',
        'order',
        'family',
        'genus',
        'species',
        'strain'
      )) %>% 
    nest(tax_lineage = c(
      'no rank',
      'superkingdom',
      'clade',
      'phylum',
      'class',
      'order',
      'family',
      'genus',
      'species',
      'strain'
    )) %>% 
    select(-xml)
}

