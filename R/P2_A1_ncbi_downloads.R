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
    # dfs <- dfs |>
    #   mutate(
    #     nuc_ss = map(get, fasta_to_DNAss),
    #     nuc_name = map(nuc_ss, names),
    #     nuc_seq = map(nuc_ss, paste)
    # )
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

# return taxonomic records for ids
# input dataframe with column called tax_id,
fetch_taxonomy <- function(data, id, chunk_size = 50){
  
  message('input data')
  # only keep unique identifiers
  df <- data |> 
    transmute(id = {{id}}) |> 
    distinct() 
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
        entrez_post(db = 'taxonomy', id =.x)
      }))
  
  message('Fetching summaries')
  pb <- progress_bar$new(total = nrow(dfs))
  dfs <- dfs |> 
    mutate(
      xml = map(token, ~{
        Sys.sleep(0.2)
        pb$tick()
        entrez_fetch(db = 'taxonomy', web_history = .x)
        })
    )
  toc()
  return(dfs)
}



# return taxonomic records for ids
# input dataframe with column called tax_id,
parse_taxonomy <- function(data, xml){
#   
  data |>
    # TODO check that map_chr works here
    mutate(
      xml =
        map_chr({{id}}, ~entrez_fetch(db = 'taxonomy', id = .x,
                                      rettype="xml")),
      tax_name =
        str_extract(xml, '<ScientificName>.*?</ScientificName>') |>
        str_remove_all('<ScientificName>|</ScientificName>'),
      tax_parent =
        str_extract(xml, '<ParentTaxId>[0-9]+') |>
        str_remove_all('<ParentTaxId>'),
      tax_lineage =
        str_extract(xml, '<Lineage>.*?</Lineage>\n') |>
        str_remove_all('<Lineage>|</Lineage>\n'),
      tax_genetic_code =
        str_extract(xml, '<GeneticCode>\n        <GCId>[0-9]+') |>
        str_remove_all('<GeneticCode>\n        <GCId>'),
    )  |>
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
      )) |>
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
    )) |>
    select(-xml)
}


# MAIN ------------------------------------------------------

## Get all cdd and prot ids first --------

# CDD search for list of RIT terms
# terms <- tibble(term = map_chr(c('A', 'B', 'C'), ~glue('INT_Rit{.x}_C_like')))
# cdd_searches <- terms |> cdd_search(term = term)
# cdd_summary <- cdd_searches |> unnest_wider(cdd_summary)
# write_rds(cdd_summary, './data/CDD/cdd_summary.rds')

# ## Link cdd records to proteins and nucleotides
# cdd_prot_ids <-
#   cdd_searches |>
#   select(-cdd_summary) |>
#   link_cdd_protein(id = cdd_id) |> 
#   unnest(prot_id)
# beep()
# write_rds(cdd_prot_ids, './data/CDD/cdd_prot_ids.rds')
cdd_prot_ids <- read_rds('./data/CDD/cdd_prot_ids.rds')


## link prot to nuc -----
# 
# ## DIDN'T WORK
# prot_ids <- cdd_prot_ids |> select(prot_id) |> distinct()
# 
# n <- nrow(cdd_prot_ids)
# chunk_size <- 100
# pb <- progress_bar$new(
#   format = "  downloading [:bar] :percent eta: :eta",
#   total = n%/%100, clear = FALSE, width= 60
#   )
# 
# cdd_prot_nuc_ids <-
#   # split prot ids into in chunks
#   tibble(
#     df = map(
#       .x = seq(1, n, chunk_size),
#       .f = ~prot_ids |> slice(.x:min(.x + chunk_size, n)))
#   ) |> 
#   # get nuc ids
#   mutate(nuc_id = map(df, ~{
#     pb$tick()
#     link_protein_nuccore(.x, id = prot_id)
#   }
#   ))

## DIDNT WORK 
# cdd_prot_nuc_ids1 <- 
#   cdd_prot_ids |>
#   link_protein_nuccore(id = prot_id) |> 
#   unnest(nuc_id)
# beep()

# cdd_prot_nuc_tax_ids <- 
#   cdd_prot_nuc_ids |> 
#   link_nuccore_taxonomy(id = nuc_id)

# write_rds(cdd_prot_nuc_tax_ids, './data/CDD/dataset_ids.rds')

# id_data <- cdd_prot_nuc_tax_ids

## Download data -------
## check number of records before committing to download

id_data <- read_rds('./data/CDD/prot_and_nuc_ids.rds')

# # get protein data and summaries
# prot_data <- id_data |>
#   fetch_data(id = prot_id, db = 'protein', chunk_size = 50) |> 
#   unnest(cols = c(id, prot_name, prot_seq)) |> 
#   transmute(prot_id = id, prot_name, prot_seq)
# beep()
# write_rds(prot_data, './data/CDD/prot_data.rds')
# 
# prot_summary <- id_data |> 
#   fetch_summaries(id = prot_id, db = 'protein', chunk_size = 50)
# 
# prot_summary |> select(-token) |> unnest_wider(summary)
# beep()
# write_rds(prot_summary, './data/CDD/prot_summary.rds')
# rm(prot_data, prot_summary, prot_ids)


##
nuc_ids <- id_data |>
  unnest(nuc_id) |> 
  select(nuc_id) |> 
  distinct()
rm(id_data)

dir.create('./data/CDD/nuc_data')


### better way to create files
fetch_nuc_datasets <- function(data, set_number) {
  message(glue('set number: {set_number}'))
  data |> 
    fetch_data(id = nuc_id, db = 'nuccore', chunk_size = 100) |> 
    parse_nuc_data() |> 
    write_rds(nuc_data12, glue('./data/CDD/nuc_data/{set_number}.rds'))
  return('./data/CDD/nuc_data/{set_number}.rds')
}

n <- nrow(nuc_ids)
nuc_sets <- 
  tibble(
    id_set = map(seq(1, n, 1000), ~nuc_ids[.x: min(.x + 1000, n),])
    ) |> 
  mutate(set_number = row_number(),
         filepath = map2(
           .x = id_set, 
           .y = set_number,
           .f = ~fetch_nuc_datasets(data = .x, set_number = .y)
           )
         )


# 
# ## for 16 individual sets....
# nuc_data1 <- nuc_ids |> slice(1:1000) |> 
#   fetch_data(id = nuc_id, db = 'nuccore', chunk_size = 100) |> 
#   parse_nuc_data()
# write_rds(nuc_data1, './data/CDD/nuc_data/1.rds', compress = 'gz')
# beep()
# rm(nuc_data1)
# 
# nuc_data2 <- nuc_ids |> slice(1001:2000) |> 
#   fetch_data(id = nuc_id, db = 'nuccore', chunk_size = 100) |> 
#   parse_nuc_data()
# write_rds(nuc_data2, './data/CDD/nuc_data/2.rds')
# beep()
# rm(nuc_data2)
# 
# nuc_data3 <- nuc_ids |> slice(2001:3000) |> 
#   fetch_data(id = nuc_id, db = 'nuccore', chunk_size = 100) |> 
#   parse_nuc_data()
# write_rds(nuc_data3, './data/CDD/nuc_data/3.rds')
# beep()
# rm(nuc_data3)
# 
# nuc_data4 <- nuc_ids |> slice(3001:4000) |> 
#   fetch_data(id = nuc_id, db = 'nuccore', chunk_size = 100)  |> 
#   parse_nuc_data()
# write_rds(nuc_data4, './data/CDD/nuc_data/4.rds')
# beep()
# rm(nuc_data4)
# 
# nuc_data5 <- nuc_ids |> slice(4001:5000) |> 
#   fetch_data(id = nuc_id, db = 'nuccore', chunk_size = 100) |> 
#   parse_nuc_data()
# write_rds(nuc_data5, './data/CDD/nuc_data/5.rds')
# beep()
# rm(nuc_data5)
# 
# nuc_data6 <- nuc_ids |> slice(5001:6000) |> 
#   fetch_data(id = nuc_id, db = 'nuccore', chunk_size = 100) |> 
#   parse_nuc_data()
# write_rds(nuc_data6, './data/CDD/nuc_data/6.rds')
# beep()
# rm(nuc_data6)
# 
# nuc_data7 <- nuc_ids |> slice(6001:7000) |> 
#   fetch_data(id = nuc_id, db = 'nuccore', chunk_size = 100) |> 
#   parse_nuc_data()
# write_rds(nuc_data7, './data/CDD/nuc_data/7.rds')
# rm(nuc_data7)
# beep()
# 
# nuc_data8 <- nuc_ids |> slice(7001:8000) |> 
#   fetch_data(id = nuc_id, db = 'nuccore', chunk_size = 100) |> 
#   parse_nuc_data()
# write_rds(nuc_data8, './data/CDD/nuc_data/8.rds')
# beep()
# rm(nuc_data8)
# 
# nuc_data9 <- nuc_ids |> slice(8001:9000) |> 
#   fetch_data(id = nuc_id, db = 'nuccore', chunk_size = 100) |> 
#   parse_nuc_data()
# write_rds(nuc_data9, './data/CDD/nuc_data/9.rds')
# rm(nuc_data9)
# beep()
# 
# nuc_data10 <- nuc_ids |> slice(9001:10000) |> 
#   fetch_data(id = nuc_id, db = 'nuccore', chunk_size = 100) |> 
#   parse_nuc_data()
# write_rds(nuc_data10, './data/CDD/nuc_data/10.rds')
# rm(nuc_data10)
# beep()
# 
# nuc_data11 <- nuc_ids |> slice(10001:11000) |> 
#   fetch_data(id = nuc_id, db = 'nuccore', chunk_size = 100) |> 
#   parse_nuc_data()
# write_rds(nuc_data11, './data/CDD/nuc_data/11.rds')
# rm(nuc_data11)
# beep()
# 
# nuc_data12 <- nuc_ids |> slice(11001:nrow(nuc_ids)) |> 
#   fetch_data(id = nuc_id, db = 'nuccore', chunk_size = 100) |> 
#   parse_nuc_data()
# beep()
# write_rds(nuc_data12, './data/CDD/nuc_data/12.rds')
# rm(nuc_data12)
