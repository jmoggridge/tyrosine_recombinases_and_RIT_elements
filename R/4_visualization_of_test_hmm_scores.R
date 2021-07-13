library(tidyverse)  
library(cluster)    
library(factoextra) 
library(umap)
library(pals)
library(patchwork)

test_df <- 
  read_rds('./results/final_model_07-05/final_prep.rds') |> 
  pull(test_prep) |> 
  pluck(1) |> 
  select(subfamily, acc, Arch1:Xer)

write_csv(test_df |> janitor::clean_names(),
          './results/test_hmm_scores.csv')

df <- test_df |> select(-subfamily, -acc) |> as.data.frame()
rownames(df) <- test_df$acc

set.seed(123)
clusters <- kmeans(df, centers = 21, nstart = 40)
str(clusters)

library(gt)
test_df |> 
  bind_cols(cluster = clusters$cluster) |> 
  count(cluster, subfamily) |> 
  arrange(cluster, desc(n)) |> 
  # mutate(cluster = paste('Cluster', cluster)) |> 
  gt::gt(auto_align = T) |> 
  tab_options(
    row_group.padding = 1L, 
    table.font.size = 11,
    data_row.padding = 1L
    ) |> 
  tab_header(
    title = "K-means clustering of HMM scores of test sequences",
    subtitle = "k = 21"
    ) |> 
  cols_align(align='left', columns = subfamily) |> 
  data_color(
    columns = cluster,
    alpha = 0.5,
    colors = scales::col_factor(
      palette = as.character(pals::alphabet2(21)),
      domain = NULL
    ))



fviz_cluster(
  clusters, 
  data = df,
  geom = c("point"),
  ellipse.alpha = 0.01,
  shape = 1,
  pointsize = 0.25,
  outlier.labelsize = 0,
  ggtheme = theme_bw(),
)



## JUST PCA -----

pc <- prcomp(df)
pcx <- pc$x |> 
  as.data.frame() |> 
  rownames_to_column(var = 'acc') |> 
  select(acc:PC3) |> 
  as_tibble() |> 
  right_join(test_df |> select(acc, subfamily))

p1 <- pcx |> 
  ggplot(aes(PC1, PC2, color = subfamily)) + 
  geom_point(shape = 1, size = 0.5, alpha = 0.4) +
  scale_color_discrete(type =  as.character(pals::alphabet2(21))) +
  guides(colour = guide_legend(
    override.aes = list(size=3, shape = 15, alpha = 1))
  ) +
  theme_bw()
p2 <- pcx |> 
  ggplot(aes(PC1, PC3, color = subfamily)) + 
  geom_point(shape = 1, size = 0.5, alpha = 0.4) +
  scale_color_discrete(type =  as.character(pals::alphabet2(21))) +
  guides(colour = guide_legend(
    override.aes = list(size=3, shape = 15, alpha = 1))
  ) +
  theme_bw()
p3 <- pcx |> 
  ggplot(aes(PC2, PC3, color = subfamily)) + 
  geom_point(shape = 1, size = 0.5, alpha = 0.4) +
  scale_color_discrete(type =  as.character(pals::alphabet2(21))) +
  guides(colour = guide_legend(
    override.aes = list(size=3, shape = 15, alpha = 1))
  ) +
  theme_bw()

p1 + p2 + p3 + guide_area() +
  plot_layout(guides = 'collect') &
  plot_annotation(title = 'PCA of HMM scores of test sequences')



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
  
  







