library(tidyverse)
library(rentrez)
library(progress)
library(tictoc)

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
cdd_link_protein <- function(data, id) {
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
protein_link_nuccore <- function(data, id){
  data |>
    unnest({{id}}) |>
    # can get db elinks from proteins to  nuccore records
    mutate(link = map(
      .x = {{id}},
      .f = ~ entrez_link(
        id = .x,
        dbfrom = "protein",
        db = "nuccore"
      )
    ),
    nuc_id = map(link, ~.x$links$protein_nuccore)) |>
    # drop any protein records without a linked nucleotide
    mutate(drop = map_lgl(nuc_id, ~{!is.null(.x)})) |> 
    filter(drop == T) |> 
    select(-c(drop, link))
}

# get linked taxonomy records for each nucleotide
# input tibble with a column of ids
nuccore_link_taxonomy <- function(data, id){
  taxonomy_ids <- data |>
    select({{id}}) |>
    distinct() |>
    mutate(
      tax_link = map(
        .x = {{id}},
        .f = ~entrez_link(
          id = .x,
          dbfrom = 'nuccore',
          db = 'taxonomy',
        )
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
    dfs <- dfs |>
      mutate(
        nuc_ss = map(get, fasta_to_DNAss),
        nuc_name = map(nuc_ss, names),
        nuc_seq = map(nuc_ss, paste)
      )
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
        entrez_fetch(db = 'taxonomy', web_history = .x)  |> 
          enframe() |> 
          unnest_wider(value)})
    )
  toc()
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

## CDD search for list of RIT terms
terms <- tibble(term = map_chr(c('A', 'B', 'C'), ~glue('INT_Rit{.x}_C_like')))
cdd_searches <- terms |> cdd_search(term = term)
cdd_summary <- cdd_searches |> unnest_wider(cdd_summary)
write_rds(cdd_summary, './data/CDD/cdd_summary.rds')

# ## Link cdd records to proteins and nucleotides
cdd_prot_ids <-
  cdd_searches |>
  select(-cdd_summary) |>
  cdd_link_protein(id = cdd_id) |> 
  unnest(prot_id)

cdd_prot_nuc_ids <- 
  cdd_prot_ids |>
  protein_link_nuccore(id = prot_id) |> 
  unnest(nuc_id)
beep()

cdd_prot_nuc_tax_ids <- 
  cdd_prot_nuc_ids |> 
  nuccore_link_taxonomy(id = nuc_id)

write_rds(cdd_prot_nuc_tax_ids, './data/CDD/dataset_ids.rds')

id_data <- 
## check number of records before commiting to download

# id_data <- 
#   read_rds('./data/CDD/prot_and_nuc_ids.rds') |> 
#   unnest(nuc_id)

prot_data <- id_data |>
  fetch_data(id = prot_id, db = 'protein', chunk_size = 50) |> 
  unnest(cols = c(id, prot_name, prot_seq)) |> 
  transmute(prot_id = id, prot_name, prot_seq)

prot_summary <- id_data |> 
  fetch_summaries(id = nuc_id, db = 'protein', chunk_size = 50)

nuc_data <- id_data |>
  fetch_data(id = nuc_id, db = 'nuccore', chunk_size = 10) |> 
  unnest(cols = c(id, nuc_name, nuc_seq)) |> 
  transmute(nuc_id = id, nuc_name, nuc_seq)

nuc_summary <- id_data |>
  fetch_summaries(id = nuc_id, db = 'nuccore', chunk_size = 50)



beepr::beep()
glimpse(try_this)


try_this |> 
  select(id, contains('name'), contains('seq')) |> 
  unnest(cols = c(id, nuc_name, nuc_seq))

id_data |> 
  transmute(id = tax_id)


id_data |> pull(tax_id) |> unique() |> length()
# 
# id_data <- 
#   read_rds('./data/CDD/prot_and_nuc_ids.rds') |> 
#   unnest(nuc_id)
# 
# nuc_df <- 
#   id_data |> 
#   select(nuc_id) |> 
#   distinct() |> 
#   slice(1:200)
# 
# nuc_dfs <- tibble(
#   id_set = map(.x = seq(1, nrow(nuc_df), 50),
#                .f = ~nuc_df |> slice(.x:min(.x+49, nrow(nuc_df))))
# ) |> 
#   mutate()
# 
# nuc_dfs <- nuc_dfs |> 
#   mutate(token = map(id_set, ~{
#     Sys.sleep(0.2)
#     entrez_post(db = 'nuccore', id =.x$nuc_id)
#     }))
# 
# nuc_data <- nuc_dfs |> 
#   mutate(get = map(token, ~{
#     entrez_fetch(db = 'nuccore', web_history = .x, rettype = 'fasta')
#   }))
# 
# nuc_dat2 <- nuc_data |> 
#   mutate(nuc_summary = map(token, ~{
#     Sys.sleep(0.2)
#     entrez_summary(db = 'nuccore', web_history = .x)  |> 
#       enframe() |> 
#       unnest_wider(value)
#   }))
# nuc_dat2
# nuc_dat2 |> 
#   unnest()
# 
# entrez_summary(db = 'nuccore', web_history = nuc_data$token[[1]]) |> 
#   enframe() |> 
#   unnest_wider(value)
# 
# # fasta to DNA set
# fasta_to_DNAss <- function(x){
#   temp <- tempfile()
#   write_file(x, temp)
#   Biostrings::readDNAStringSet(temp)
# }
# 
# if (db == 'nuccore'){
#   
# }
# nuc_data |> 
#   mutate(nuc_ss = map(get, fasta_to_DNAss),
#          nuc_name = map(nuc_ss, names),
#          nuc_seq = map(nuc_ss, paste)
#          )
# 
# nuc_data |> unnest_wider(id_set)