# P2_AXX Visualizations

rits <- read_rds('results/non_redundant_rits.rds')

# From 1087 nucleotide records,
# hard to figure out what is partial or complete / contigs v scaffolds
# finished vs draft.
# - ~160 complete genome/plasmid
# - 38 records indicate draft genome
# rest are contig-level sequences.
rits |> 
  select(nuc_id, completeness, sourcedb, genome, nuc_name) |> 
  distinct() |> 
  mutate(
    st = map_chr(
      nuc_name,
      ~{
        tolower(.x) |> 
          str_extract_all(
            'draft|complete genome|wgs|contig|scaffold|complete plasmid'
          ) |> 
          unique() |> 
          paste0(collapse = ' ')
      }),
    st2 = glue('{st}_{completeness}_{genome}_{sourcedb}')
  ) |> 
  # View()
  count(st2) |> 
  print(n=30)


rits |> 
  select(nuc_id, completeness, sourcedb, slen, genome, nuc_name) |> 
  distinct() |> 
  pull(slen) |> 
  summary()

# can see that 'chromosome' mostly genomes but there are far more contigs ~3kb
rits |> 
  select(nuc_id, completeness, sourcedb, slen, genome, nuc_name) |> 
  distinct() |> 
  ggplot(aes(slen, fill = genome)) +
  geom_histogram() +
  scale_x_log10() +
  labs(
    fill = 'ncbi label',
    title = 'how long are the contigs or scaffolds scaffolds?')






















