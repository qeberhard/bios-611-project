library(ggplot2)
library(stringr)
library(readr)
library(dplyr)

#Here we are going to attempt to cluster the RNAs based on their k-mer similarity
transcriptome1 <- read_csv("source_data/5mers_expressed_txpts_1.csv.gz")
transcriptome2 <- read_csv("source_data/5mers_expressed_txpts_2.csv.gz")
transcriptome <- rbind(transcriptome1, transcriptome2)

transcriptome_wide <- transcriptome %>% 
  arrange()
arrange(time) %>%
  pivot_wider(id_cols=c("trial","label"),
              names_from="time",
              values_from="V") %>%
  mutate(index=1:nrow(.)) %>%
  arrange(runif(nrow(.)));
#ggplot(transcriptome, aes(ACACA,GUGUG)) + geom_point() + coord_equal()

results <- kmeans(transcriptome_wide, centers=3);
data_clusters <- transcriptome %>% mutate(cluster=results$cluster)
png("figures/kmeans_cluster_3_centered.png")
ggplot(data_clusters, aes(x,y)) + geom_point(color=factor(data_clusters$cluster)) + coord_equal()
dev.off()


results <- kmeans(transcriptome_wide, centers=4);
data_clusters <- transcriptome %>% mutate(cluster=results$cluster)
png("figures/kmeans_cluster_4_centered.png")
ggplot(data_clusters, aes(x,y)) + geom_point(color=factor(data_clusters$cluster)) + coord_equal()
dev.off()

results <- kmeans(transcriptome_wide, centers=5);
data_clusters <- transcriptome %>% mutate(cluster=results$cluster)
png("figures/kmeans_cluster_5_centered.png")
ggplot(data_clusters, aes(x,y)) + geom_point(color=factor(data_clusters$cluster)) + coord_equal()
dev.off()

results <- kmeans(transcriptome_wide, centers=6);
data_clusters <- transcriptome %>% mutate(cluster=results$cluster)
png("figures/kmeans_cluster_6_centered.png")
ggplot(data_clusters, aes(x,y)) + geom_point(color=factor(data_clusters$cluster)) + coord_equal()
dev.off()