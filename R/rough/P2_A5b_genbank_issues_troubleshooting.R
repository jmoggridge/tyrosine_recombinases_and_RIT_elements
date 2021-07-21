# these nucleotides probably contain RIT element bc have 3 CDD proteins.
# TODO go back to P2_A4 and find out why duplicated rows here
three_ints <- read_rds('./data/CDD/ids_w_three_integrases.rds') |> 
  arrange(slen) |> 
  distinct()
glimpse(three_ints)

# create an index for nuc ids to their genbank records
genbank_files_index <- 
  tibble(file = Sys.glob('./data/CDD/RIT_gbk_[0-9].rds')) |> 
  mutate(nuc_id = map(file, ~ read_rds(.x) |> pull(nuc_id))) |> 
  unnest(nuc_id) |> 
  group_by(nuc_id) |> 
  filter(row_number() == 1) |> 
  ungroup()


have_genbank <- three_ints |> 
  filter(nuc_id %in% genbank_files_index$nuc_id)

# these are ids that genbank records aren't present for
no_genbank <- 
  anti_join(three_ints, have_genbank) |> pull(nuc_id) |> unique()
no_genbank

## retry these 
missing <- read_rds('./data/CDD/rit_finder/missing_results.rds')

# missing_genbank |> 
#   map2(file, ~open_genba())....
