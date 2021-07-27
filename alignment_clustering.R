
library(Biostrings)
library(DECIPHER)
library(tidyverse)
library(ape)

# filtered rits data from P2-A7;
nr_rits <- 
  read_rds('results/non_redundant_rits.rds') |> 
  mutate(
    rit_dna = str_to_upper(rit_dna),
    clade = ifelse(clade == 'Bacteria incertae sedis', 'unclassified Bacteria', clade)
  ) 
glimpse(nr_rits)  

# set of distinct sequences
uniq_rit_dna <- nr_rits |> 
  select(rit_dna, rit_id) |> 
  distinct()

# any sequences unusually long or short? 
# the Nitrospirae sequence seems short and one in the PVC group
# bimodal @ ~3100 & ~3500
nr_rits |> 
  ggplot(aes(rit_length)) + 
  geom_histogram(show.legend = F) +
  facet_wrap(~phylum, scales = 'free_y')


# investigate base composition:
# need to fix non-IUPAC bases? no.
str_split('ACGTNUDBSYRWKMV-', '', simplify = T) |>
  c('[^ACGTNUDBSYRWKMV\\-]') |> 
  set_names() |> 
  map(~str_count(uniq_rit_dna$rit_dna, .x) / nchar(uniq_rit_dna$rit_dna)) |> 
  enframe(name = 'base', value = 'prop') |> 
  unnest(prop) |> 
  group_by(base) |> 
  summarize(across(.cols = c(prop), 
                   .fns = list(min = min, mean = mean, med = median, mx = max)))

# # A tibble: 17 x 5
# base                    prop_min   prop_mean prop_med  prop_mx
# <chr>                      <dbl>       <dbl>    <dbl>    <dbl>
#   1 "-"                        0     0              0     0       
# 2 "[^ACGTNUDBSYRWKMV\\-]"    0     0              0     0       
# 3 "A"                        0.127 0.237          0.204 0.394   
# 4 "B"                        0     0              0     0       
# 5 "C"                        0.128 0.270          0.297 0.388   
# 6 "D"                        0     0              0     0       
# 7 "G"                        0.153 0.263          0.282 0.348   
# 8 "K"                        0     0.000000609    0     0.000558
# 9 "M"                        0     0              0     0       
# 10 "N"                        0     0.0000268      0     0.0214  
# 11 "R"                        0     0              0     0       
# 12 "S"                        0     0              0     0       
# 13 "T"                        0.147 0.229          0.218 0.330   
# 14 "U"                        0     0              0     0       
# 15 "V"                        0     0              0     0       
# 16 "W"                        0     0.000000304    0     0.000279
# 17 "Y"                        0     0              0     0       




# Alignment and distance for unique RIT seqs ------

dna_ss <- uniq_rit_dna$rit_dna |> 
  set_names(nm = uniq_rit_dna$rit_id) |> 
  Biostrings::DNAStringSet(use.names = T)
dna_ss
 
Biostrings::writeXStringSet(dna_ss, 'data/rit_dna.fasta')

dna_aligned <- DECIPHER::AlignSeqs(
  myXStringSet = dna_ss, 
  verbose = T, 
  iterations = 5)
beepr::beep()
Biostrings::writeXStringSet(dna_aligned, 'data/rits_alignged.fasta')

# DECIPHER::BrowseSeqs(dna_aligned)
dna_bin <- ape::as.DNAbin(dna_aligned)

# get 'raw' distance matrix
dist_df <- 
  ape::dist.dna(x = dna_bin, 
                model = 'raw', 
                pairwise.deletion = T, 
                as.matrix = T) |> 
  as_tibble() |> 
  mutate(rit_id = names(dna_bin)) |> 
  pivot_longer(cols = -rit_id, values_to = 'distance') |> 
  filter(rit_id != name)

# see furthest pairs ~ 0.7
dist_df |> 
  arrange(desc(distance)) |> 
  print(n = 200)

# see closest pairs; many 0's
dist_df |> 
  arrange(distance) |> 
  print(n = 200)

## Ones that are 0 distance can be just a subsequence - misleading.
# reduced_groups <- dist_df |> 
#   filter(distance == 0) |> 
#   rowwise() |> 
#   mutate(rit1 = min(rit_id, name),
#          rit2 = max(name, rit_id)) |> 
#   select(rit1, rit2) |> 
#   distinct() |> 
#   group_by(rit1) |> 
#   summarize(same_group = list(rit2)) |> 
#   mutate(n = map_int(same_group, length)) |> 
#   arrange(desc(n)) |> unnest(same_group)
# 
# reduced_groups[which(reduced_groups$same_group %in% reduced_groups$rit1),]

# get 'raw' distance
dist_n_df <- 
  ape::dist.dna(x = dna_bin, model = 'N') |> 
  as.matrix() |> 
  as_tibble() |> 
  mutate(rit_id = names(dna_bin)) |> 
  pivot_longer(cols = -rit_id, values_to = 'distance') |> 
  rowwise() |>
  mutate(rit1 = min(rit_id, name),
         rit2 = max(name, rit_id)) |>
  select(rit1, rit2, distance) |>
  distinct()
dist_n_df |> 
  arrange(desc(distance)) |> 
  print(n = 200)
dist_n_df |> 
  arrange(distance) |> 
  print(n = 200)

distance_matrix <- DECIPHER::DistanceMatrix(dna_aligned, type='dist', includeTerminalGaps = T)
class(distance_matrix)

# Cluster seqs ------

clusters_upgma <- DECIPHER::IdClusters(
  myXStringSet = dna_aligned,
  myDistMatrix = distance_matrix,
  method = 'UPGMA',  
  cutoff = c(0.3, 0.15, 0.1, 0.05, 0.01,  0.005, 0.001), 
  type = "clusters", 
  showPlot = F, 
  verbose = TRUE
  )
clusters_upgma

clusters_nj <- DECIPHER::IdClusters(
  myXStringSet = dna_aligned,
  myDistMatrix = distance_matrix,
  method = 'NJ',  
  cutoff = c(0.3, 0.15, 0.1, 0.05, 0.01,  0.005, 0.001), 
  type = "clusters", 
  showPlot = F, 
  verbose = TRUE
) 
clusters_nj


# NJ tree
nj <- clusters <- DECIPHER::IdClusters(
  myXStringSet = dna_aligned,
  myDistMatrix = distance_matrix,
  method = 'NJ',  
  type = "dendrogram", 
  verbose = TRUE
)
# upgma_tree 
upgma <- DECIPHER::IdClusters(
  myXStringSet = dna_aligned,
  myDistMatrix = distance_matrix,
  method = 'UPGMA',  
  type = "dendrogram", 
  verbose = TRUE
)
par(cex=1, mar=c(5, 8, 4, 1))
plot(upgma,
     leaflab = 'none',
     main = 'UPGMA tree of RIT coding region dna sequences',
     horiz = T)
par(cex=1, mar=c(5, 8, 4, 1))
plot(nj,
     leaflab = 'none',
     main = 'NJ tree of RIT coding region dna sequences',
     horiz = T)



## cluster based on protein sequences...

# this rit is unusually long! and aligns relatively poorly 
# nr_rits |> 
#   filter(rit_id == 'RIT_171') |> 
#   View()


rits_clustered <- 
  tibble(clusters_upgma) |> 
  mutate(rit_id = rownames(clusters_upgma)) |> 
  right_join(nr_rits, by = 'rit_id')

rits_clustered |> 
  group_by(cluster0_005UPGMA) |> 
  count(nuc_id) |> 
  arrange(desc(n))


nr_rits |> filter(rit_id %in% c('RIT_208', 'RIT_709')) 

tibble(rits_clustered)
rits_clustered |> count(cluster0_3UPGMA) |> nrow()
clusters_upgma |> count(cluster0_15UPGMA) |> nrow()
clusters_upgma |> count(cluster0_1UPGMA) |> nrow()
clusters_upgma |> count(cluster0_05UPGMA) |> nrow()
clusters_upgma |> count(cluster0_01UPGMA) |> nrow()
clusters_upgma |> count(cluster0_005UPGMA) |> nrow()
clusters_upgma |> count(cluster0_001UPGMA) |> nrow()

clusters_nj |> count(cluster0_3NJ) |> nrow()
clusters_nj |> count(cluster0_15NJ) |> nrow()
clusters_nj |> count(cluster0_1NJ) |> nrow()
clusters_nj |> count(cluster0_05NJ) |> nrow()
clusters_nj |> count(cluster0_01NJ) |> nrow()
clusters_nj |> count(cluster0_005NJ) |> nrow()
clusters_nj |> count(cluster0_001NJ) |> nrow()


length(dna_ss)




## Active elements for terminal repeats -----

# select nucleotides with > 1 copy; easiest and will have enough diversity of examples?

nr_rits |> 
  group_by(nuc_id, rit_id) |>  # identical sequences 
  mutate(rit_count = length(rit_dna)) |> 
  ungroup() |> 
  count(rit_count) |> 
  arrange(desc(n))

multiple_rit_copies <- nr_rits |> 
  group_by(nuc_id, rit_id) |>  # identical coding sequences could have different upstream and downstream regions?
  mutate(rit_count = length(rit_dna)) |> 
  filter(rit_count > 1) |> 
  select(rit_count, rit_id, nuc_id, rit_dna_upstream, rit_dna_downstream) |> 
  distinct() |> 
  mutate(
    head_tail = glue('{rit_dna_upstream}{revcomp(rit_dna_downstream)}'),
    head_tail_name = glue('{rit_id}_{nuc_id}')
    )
  


head_tail <- DNAStringSet(multiple_rit_copies$head_tail)
names(head_tail) <- multiple_rit_copies$head_tail_name
DECIPHER::DetectRepeats(head_tail, type = 'both')

