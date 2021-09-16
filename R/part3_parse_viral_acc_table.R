# Get virus genomes using genome tables from
# https://www.ncbi.nlm.nih.gov/genomes/GenomesGroup.cgi?taxid=10239&host=bacteria
# and 
# https://www.ncbi.nlm.nih.gov/genomes/GenomesGroup.cgi?taxid=10239&host=archaea


library(tidyverse)
library(janitor)
library(lubridate)
library(rentrez)
library(progress)

# parse tables for bacterial and archaeal viruses
all_entries <- 
  bind_rows(
    read_tsv('data/ncbi_viral_genomes/ncbi_bacteria_viruses.tbl', skip = 1),
    read_tsv('data/ncbi_viral_genomes/ncbi_archaea_viruses.tbl', skip = 1)
  ) |> 
  clean_names() |> 
  mutate(genome_length = str_extract(genome_length, '[0-9]+') |> as.numeric(),
         genome_neighbors = str_extract(genome_length, '[0-9]+') |> as.numeric(),
         date_completed = mdy(date_completed),
         date_updated = mdy(date_updated),
         number_of_segments =  as.numeric(number_of_segments)
  ) |> 
  mutate(across(where(is.character), as_factor)) |> 
  distinct()

glimpse(all_entries)


# genomes with only a single accession
single_segments <- all_entries |> 
  filter(!is.na(accession)) |> 
  filter(accession != '-')


# some genomes have multiple segments under different accessions
# needs to be parsed differently and joined with their header row (with no accession)
multiple_segments <- all_entries |> 
  filter(is.na(accession)) |> 
  select(genome) |> 
  separate(genome, sep = '   ',
           into = c('genome', 'genome_length', 'accession',
                    'number_of_proteins', 'genome_neighbors')
           ) |> 
  mutate(accession = str_trim(accession)) |> 
  mutate(across(
    .cols = c(genome_length, number_of_proteins, genome_neighbors),
    .fns = ~str_extract(.x, '[0-9]+') |> as.numeric()
    )
  ) 


# keep the metadata from the ones with multiple accessions
multiple_segments <- 
  all_entries |> 
  filter(accession == '-') |> 
  select(-c(accession, genome_length, number_of_proteins, genome_neighbors)) |> 
  full_join(
    multiple_segments,
    by = 'genome'
  )
glimpse(multiple_segments)


# join the single and multiple segment genomes
virus_tbl <- 
  full_join(
    single_segments,
    multiple_segments
  ) 

glimpse(virus_tbl)  
summary(virus_tbl)

sum(is.na(virus_tbl$accession))
sum(virus_tbl$accession == '-')

rm(all_entries, single_segments, multiple_segments)




fetch_viral_genome <- function(accession){
  search <- rentrez::entrez_search(db =  'nuccore', term = accession, use_history = T)
  genbank_xml <- rentrez::entrez_fetch(
    db = 'nuccore', 
    web_history = search$web_history, 
    retmode = 'xml',
    rettype = 'gb')
  return(genbank_xml)
}


acc <- virus_tbl |> 
  select(genome, accession) |> 
  slice(1:5) 

pb <- make_pb(nrow(acc))

acc |> 
  mutate(xml = map(accession, ~{
    xml <- fetch_viral_genome(.x)
    pb$tick()
    return(xml)
    }
    ))



r <- parse_genbank(rcd)
glimpse(r)
r





