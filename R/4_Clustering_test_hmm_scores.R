library(tidyverse)  
library(cluster)    
library(factoextra) 
library(umap)
library(pals)
library(patchwork)
library(beepr)
library(clustertend)
test_df <- 
  read_rds('./results/final_model_07-11/prepped_data.rds') |> 
  pull(test_prep) |> 
  pluck(1) |> 
  select(subfamily, acc, Arch1:Xer) 
test_df |> count(subfamily) |> print(n=25)

df <- test_df |> 
  select(-subfamily, -acc) |>
  as.data.frame()

rownames(df) <- test_df$acc

# write_csv(test_df |> janitor::clean_names(),
#           './results/test_hmm_scores.csv')

# hopkins statistic: is the data normally distributed
# is H < 0.05 reject null hypothesis,
# data are clusterable: H = 0.03
set.seed(123)
# hopkins(df, n=1000)


# downsampled...
ds_df <- test_df |> 
  group_by(subfamily) |> 
  slice_sample(n=50, replace = F) |> 
  ungroup() |> 
  select(-subfamily, -acc) |>
  as.data.frame()

# factoextra::get_clust_tendency(data = df, n = 100)

# Dunn index for Ward clustering
ward_dunn_index <- NbClust::NbClust(df, distance="euclidean", 
                               min.nc=10, max.nc=40, 
                               method="ward.D2", index="dindex")

tibble(ward_dunn_index) |> 
  unnest(cols = c(ward_dunn_index)) |> 
  mutate(k = names(ward_dunn_index)) |> 
  ggplot(aes(k, ward_dunn_index)) +
  geom_path(aes(group = 1)) +
  geom_point() +
  theme_bw() +
  labs(title = 'Ward (d2) method clustering Dunn index by k')
  
beep()

# Dunn index for Kmeans 
kmeans_dunn_index <- NbClust::NbClust(df, distance="euclidean", 
                               min.nc=10, max.nc=40, 
                               method="kmeans", index="dindex")

tibble(kmeans_dunn_index) |> 
  unnest(cols = c(kmeans_dunn_index)) |> 
  mutate(k = names(kmeans_dunn_index)) |> 
  ggplot(aes(k, kmeans_dunn_index)) +
  geom_path(aes(group = 1)) +
  geom_point() +
  theme_bw() +
  labs(title = 'K-means clustering Dunn index by k')

beep()


# Ward d2 clustering method (1963)
set.seed(234)
ward_indexes <- NbClust::NbClust(
  ds_df, 
  distance="euclidean", 
  min.nc=10, 
  max.nc=40, 
  method = "ward.D2", 
  index = 'all')

beep()
# ward.d2:
# * According to the majority rule, the best number of clusters is  21 
# 
ward_indexes_df <- 
  ward_indexes$Best.nc |> t() |>
  as_tibble() |> 
  mutate(index = colnames(ward_indexes$Best.nc)) |> 
  filter(!is.na(Value_Index))

kableExtra::kable(ward_indexes_df, format = 'simple')

# Kmeans method

set.seed(234)
kmeans_indexes <- NbClust::NbClust(
  ds_df, 
  distance="euclidean", 
  min.nc=10, 
  max.nc=40, 
  method = "kmeans", 
  index = 'all')

beep()

kmeans_indexes_df <- 
  kmeans_indexes$Best.nc |> t() |>
  as_tibble() |> 
  mutate(index = colnames(kmeans_indexes$Best.nc)) |> 
  filter(!is.na(Value_Index))

kableExtra::kable(kmeans_indexes_df, format = 'simple')




# within ss
fviz_nbclust(ds_df, kmeans, method = "wss", k.max = 30) +
  theme_bw() +
  labs(title = NULL)
# silhouette
fviz_nbclust(ds_df, kmeans, method = "silhouette", k.max = 30) +
  theme_bw() +
  labs(title = NULL)
# gap statistic
fviz_nbclust(ds_df, kmeans, method = "gap_stat", k.max = 30, ) +
  theme_bw() +
  labs(title = NULL)

beep()

# make choice of k...
final <- kmeans(df, centers = 13, nstart = 50)
final2 <- kmeans(df, centers = 29, nstart = 50)

fviz_cluster(
  final,
  data = df,
  geom = c("point"),
  show.clust.cent = T,
  ellipse = F,
  ellipse.alpha = 0.00,
  shape = 16,
  pointsize = 0.5, alpha = 0.25,
  outlier.labelsize = 0,
  ggtheme = theme_classic(),
  main = 'k-means, k = 13'
)

fviz_cluster(
  final2,
  data = df,
  geom = c("point"),
  show.clust.cent = T,
  ellipse = F,
  ellipse.alpha = 0.00,
  shape = 16,
  pointsize = 0.5, alpha = 0.25,
  outlier.labelsize = 0,
  ggtheme = theme_classic(),
  main = 'k-means, k = 29'
)

# 
# set.seed(123)
# clusters <- kmeans(df, centers = 21, nstart = 40)
# str(clusters)
# 
# library(gt)
# test_df |> 
#   bind_cols(cluster = clusters$cluster) |> 
#   count(cluster, subfamily) |> 
#   arrange(cluster, desc(n)) |> 
#   # mutate(cluster = paste('Cluster', cluster)) |> 
#   gt::gt(auto_align = T) |> 
#   tab_options(
#     row_group.padding = 1L, 
#     table.font.size = 11,
#     data_row.padding = 1L
#     ) |> 
#   tab_header(
#     title = "K-means clustering of HMM scores of test sequences",
#     subtitle = "k = 21"
#     ) |> 
#   cols_align(align='left', columns = subfamily) |> 
#   data_color(
#     columns = cluster,
#     alpha = 0.5,
#     colors = scales::col_factor(
#       palette = as.character(pals::alphabet2(21)),
#       domain = NULL
#     ))


## JUST PCA -----

pc <- prcomp(df)

pcx <- pc$x |> 
  as.data.frame() |> 
  rownames_to_column(var = 'acc') |> 
  select(acc:PC3) |> 
  as_tibble() |> 
  right_join(test_df |> select(acc, subfamily))

pcx |> count(subfamily) |> print(n=25)

p1 <- pcx |> 
  ggplot(aes(PC1, PC2, color = subfamily)) + 
  geom_point(shape = 1, size = 0.2, alpha = 0.25) +
  scale_color_discrete(type =  as.character(pals::alphabet2(21))) +
  guides(colour = guide_legend(
    override.aes = list(size=3, shape = 15, alpha = 1))
  ) +
  theme_bw()
p1
p2 <- pcx |> 
  ggplot(aes(PC1, PC3, color = subfamily)) + 
  geom_point(shape = 1, size = 0.2, alpha = 0.25) +
  scale_color_discrete(type =  as.character(pals::alphabet2(21))) +
  guides(colour = guide_legend(
    override.aes = list(size=3, shape = 15, alpha = 1))
  ) +
  theme_bw()
p3 <- pcx |> 
  ggplot(aes(PC2, PC3, color = subfamily)) + 
  geom_point(shape = 1, size = 0.2, alpha = 0.25) +
  scale_color_discrete(type =  as.character(pals::alphabet2(21))) +
  guides(colour = guide_legend(
    override.aes = list(size=3, shape = 15, alpha = 1))
  ) +
  theme_bw()


p1 + p2 + p3 + plot_spacer() +
  plot_layout(guides = 'collect') &
  plot_annotation(title = 'PCA of HMM scores of test sequences')



rot_df <- pc$rotation |> 
  as_tibble() |> 
  mutate(subfamily = as_factor(rownames(pc$rotation))) |> 
  select(PC1:PC3, subfamily)

# 
# my_biplot <- function(d1, d2){
# 
#   pcx |> 
#     mutate(x = {{d1}}, y= {{d2}}) |> 
#     ggplot(aes(x, y, color = subfamily)) +
#     geom_point(shape = 1, size = 0.5, alpha = 0.15) +
#     geom_segment(
#       data = rot_df,
#       aes(x*550, y*550, xend = 0, yend= 0),
#       color = 'black'
#     ) +
#     ggrepel::geom_text_repel(
#       data = rot_df |> mutate(x={{d1}}, y = {{d2}}),
#       aes(x*550, y*550, 
#           color = subfamily, label = subfamily),
#       max.overlaps = 50
#     ) +
#     scale_color_discrete(type =  as.character(pals::alphabet2(21))) +
#     guides(colour = guide_legend(
#       override.aes = list(size=3, shape = 15, alpha = 1))
#     ) +
#     theme_classic()
# }
# 
# plts <- map2(list(1,1,2), list(2,3,3), ~my_biplot(.x, .y))
# 
# plts[[1]] + plts[[2]] + plts[[3]] & plot_layout(guides = 'collect')
# 

# GGally::ggpairs(rot_df, cardinality_threshold = 25)

## UMAP ----

test_umap <-  umap(df)
str(test_umap)

# umap_df <- 
umap_df <- test_umap$layout |> 
  as.data.frame() |> 
  rownames_to_column(var = 'acc') |> 
  as_tibble() |> 
  right_join(test_df |> select(acc, subfamily), by = "acc")

umap_df |> 
  ggplot(aes(V1, V2, color = subfamily)) +
  geom_point(size = 1, shape = 1, alpha = 0.7) + 
  theme_bw() +
  # scale_color_viridis_d() +
  scale_color_discrete(type =  as.character(pals::alphabet2(21))) +
  guides(colour = guide_legend(
    override.aes = list(size=5, shape = 15, alpha = 1))
    ) +
  labs(title = 'UMAP embedding of HMM scores for test sequences')
  
  







