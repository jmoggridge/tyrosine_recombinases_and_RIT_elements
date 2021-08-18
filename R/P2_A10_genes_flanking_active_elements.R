# Flanking genes for active elements

rits <- read_rds('results/nr_rits_clustered_IRs.rds')

# flanking genes
flanks <- rits |> 
  select(rit_id, nuc_id, tax_id, protein_df) |> 
  unnest(cols = c(protein_df)) |> 
  filter(prot_pos %in% c('lflank', 'rflank')) |> 
  filter(!is.na(prot_id))
glimpse(flanks)  

flanks |> 
  ggplot(aes(prot_pred)) +
  geom_bar() +
  coord_flip() +
  facet_wrap(~prot_pos)


# check which elements seem to be right beside each other 
has_extra_ritC <- rits |> 
  unnest(cols = c(protein_df)) |> 
  filter(prot_pos %in% c('lflank', 'rflank') & prot_pred == 'RitC') 

rits |> 
  filter(paste(rit_id, nuc_id) %in% paste(has_extra_ritC$rit_id, has_extra_ritC$nuc_id)) |> 
  unnest(protein_df) |> 
  View()





# ## what is this sequence doing here?
# read_rds('data/CDD/tax_data_fixed.rds') |> 
#   filter(tax_id =='1232866') |> 
#   unnest() |> 
#   View()
