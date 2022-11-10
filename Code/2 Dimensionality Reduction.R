################################################################################
# Setup
################################################################################
# Libraries
here::i_am("Code/2 Dimensionality Reduction.R")
library(here)       # Locate files relative to this one
library(umap)       # Calculate UMAP embeddings
library(ggplot2)    # Visualisation of UMAP embeddings / heatmaps
library(progress)   # Progress bar for creating umap plots
library(ggrepel)    # Label any outliers
library(patchwork)  # Plot umap options to pdf in a grid
library(GGally)     # Facet plot
library(tidyverse)  # Data manipulation stuff
library(clValid)    # Dunn's index (choose # clusters)
library(cluster)    # Silhouette index
library(dendextend) # Changing dendrogram aesthetics

# Read in Data
cell_df <- readRDS(here("RNA Splicing Data", "Cell Data.RDS"))
tissue_df <- readRDS(here("RNA Splicing Data", "Tissue Data.RDS"))
bio_source_all <- readRDS(here("RNA Splicing Data", "Bio_source_all.RDS"))
tissue_source <- readRDS(here("RNA Splicing Data", "Tissue source.RDS"))
################################################################################
# Functions
################################################################################
# Compare various values of n_neighbours and min)dist
create_umap_plots <- function(data, color, seed, n_components){

  # Compare choices of n and d
  n_opts <- seq(from = 5, to = 40,by = 5)
  min_dists <- c(0.001, 0.01,0.05,0.1)

  # Setup plots
  plots <- list(); plots2 <- list()
  i = 1

  # Setup progress bar
  pb <- progress_bar$new(format = "[:bar] :percent ETR: :eta",
                         total = length(n_opts)*length(min_dists),
                         complete = "-", incomplete = " ", current = ">",
                         clear = FALSE, width = 100)

  for(n in n_opts){
    for(d in min_dists){

      # Set UMAP configuration
      umap.config <- umap.defaults
      umap.config$n_neighbors <- n
      umap.config$min_dist <- d
      umap.config$n_components <- n_components

      # Compute UMAP
      set.seed(seed)
      u <- umap(data, config = umap.config, random_state = seed, transform_seed = seed)

      plot.df <- data.frame(u$layout, "col" = color)
      colnames(plot.df)

      # Create plot of results
      p1 <-  ggplot(plot.df) +
        geom_point(aes(x = X1, y = X2, col = factor(col)), show.legend = FALSE) +
        ggtitle(paste("n: ", n, ", d: ", d))
      p2 <- ggplot(plot.df) +
        geom_point(aes(x = X3, y = X1, col = factor(col)), show.legend = FALSE) +
        ggtitle(paste("n: ", n, ", d: ", d))

      # Add to list of plots
      plots[[i]] <- p1
      plots2[[i]] <- p2
      i = i + 1

      pb$tick()
    }
  }

  return(list(plots, plots2))
}

# Create inside a function for reproducibility (when compared to plots from create_umap_plots)
create_umap_single <- function(data, n_neighbors, min_dist, n_components, seed){
  # Set up configuration
  umap.config <- umap.defaults; umap.config$n_neighbors <- n_neighbors
  umap.config$min_dist <- min_dist; umap.config$n_components <- n_components

  # Run umap
  set.seed(seed)
  u <- umap(data, config = umap.config, random_state = seed, transform_seed = seed)
  return(u)
}

create_heatmap <- function(clusters, source){
  hmmat <- table(source, clusters)
  h <- heatmap(hmmat, col = hcl.colors((max(hmmat) + 1), "Reds 2", rev = TRUE))
  legend(x="topleft", legend=c(1:max(hmmat)),
         fill=hcl.colors(max(hmmat), "Reds 2", rev = TRUE))
  return(h)
}

create_ggplot_heatmap <- function(clusters, source){
  hm_df <- expand.grid(cluster = unique(clusters), source = sort(unique(source)))
  count <- as.numeric(table(clusters, source))
  hm_df$count <- ifelse(count != 0, count, NA)
  ggplot(hm_df, aes(x=cluster, y=source, fill=count)) +
    geom_tile(color="white", size = 0.25) +
    scale_fill_gradient(low="gold", high="darkorchid", na.value="white") + theme_minimal() +
    geom_text(label = ifelse(count > 0, count, ""))
}
#
#
# ################################################################################
# # Cell Data - K-means clustering
# ################################################################################
# kmeans_cell_list <- list()
# for(k in 1:30){
#   cat(k)
#   kmeans_cell_list[[k]] <- kmeans(cell_df, centers = (k+1),  nstart = 10, iter.max = 100)
# }
#
# bet_p <- c(); w_ss <- c(); tot_withinss <- c()
#
#
# for(k in 1:30) {
#   bet_p <- c(bet_p, kmeans_cell_list[[k]]$betweenss/kmeans_cell_list[[k]]$totss)
#   tot_withinss <- c(tot_withinss,  kmeans_cell_list[[k]]$tot.withinss)
#   w_ss <- c(w_ss,kmeans_cell_list[[k]]$withinss)
# }
#
# plot(2:31, tot_withinss, type = "l")
# abline(v = 4, lty = 3)
#
#
# silhouette_score <- function(kmeans, df){
#   ss <- silhouette(kmeans$cluster, dist(df))
#   mean(ss[, 3])
# }
#
# library(cluster)
# ss <- c()
# for(k in 1:30){
#   cat(k)
#   ss <- c(ss, silhouette_score(kmeans_cell_list[[k]], cell_df))
# }
#
# plot(2:31, type='b', ss, xlab='Number of clusters', ylab='Average Silhouette Scores', frame=FALSE)
#
# saveRDS(kmeans_cell_list, here("DR Files", "kmeans_cell_list.RDS"))
#
# #
# kmeans_cell_list <- readRDS(here("DR Files", "kmeans_cell_list.RDS"))
# cell_kmeans2 <- kmeans(cell_df, centers = 2,  nstart = 10, iter.max = 100)
# cell_kmeans <- kmeans(cell_df, centers = 4,  nstart = 10, iter.max = 100)
#
# cell_kmeans_df <- data.frame(cell_id = names(cell_kmeans$cluster),
#                              cluster = cell_kmeans$cluster)
# cell_kmeans_df <- left_join(cell_kmeans_df, bio_source_all)
# create_ggplot_heatmap(cell_kmeans_df$cluster, cell_kmeans_df$Biological_source)
################################################################################
# Cell - Hierarchical Clustering
################################################################################
d <- dist(cell_df)
h <- hclust(d, method = "ward.D")  # 4 clusters under complete, 2~4 under ward.D
plot(h)

ch <- cutree(h, k = 4)

h_cluster_df <- left_join(data.frame(RNA_number_id = names(ch),
                                     cluster = ch),
                          bio_source_all)
create_ggplot_heatmap(h_cluster_df$cluster, h_cluster_df$Biological_source)

################################################################################
# Cell - UMAP Reduction
################################################################################
# Create plots
# Takes ~2 minutes to run
h_plots <- create_umap_plots(cell_df, ch, 1046, 3)
#k_plots <- create_umap_plots(cell_df, cell_kmeans$cluster, 5454, 3)

plots <- h_plots[[1]]; plots2 <- h_plots[[2]]
# Compare plots
split <- length(plots)/2

pdf(here("DR Files", "Cell UMAP plots h_1.pdf"), 20, 10)
patchwork::wrap_plots(plots[1:split], nrow = 4, ncol = 4)
patchwork::wrap_plots(plots[(split + 1):length(plots)], nrow = 4, ncol = 4)
dev.off()

pdf(here("DR Files", "Cell UMAP plots h_2.pdf"), 20, 10)
patchwork::wrap_plots(plots2[1:split], nrow = 4, ncol = 4)
patchwork::wrap_plots(plots2[(split + 1):length(plots2)], nrow = 4, ncol = 4)
dev.off()

# Get umap with selected configuration
u <- create_umap_single(cell_df, 20, 0.001, 3, 1046)
ggpairs(data.frame(u$layout, "clust" = ch), columns = 1:3, aes(color = factor(clust)),
               upper = list(continuous = "blank"), diag = list(continuous = "blankDiag"))

# library(rgl)
# plot3d(x = u$layout[,1], y = u$layout[,2], z = u$layout[,3], col = ch)

cell_umap <- u$layout
saveRDS(cell_umap, here("RNA Splicing Data", "Cell UMAP.RDS"))


################################################################################
# Tissue - K-means Clustering
################################################################################
# kmeans_tissue_list <- list()
# for(k in 1:30){
#   cat(k)
#   kmeans_tissue_list[[k]] <- kmeans(tissue_df, centers = (k+1),  nstart = 10, iter.max = 100)
# }
# saveRDS(kmeans_tissue_list, "kmeans_tissue_list.RDS")
#
#
# bet_p <- c(); w_ss <- c(); tot_withinss <- c()
#
# for(k in 1:30) {
#   bet_p <- c(bet_p, kmeans_tissue_list[[k]]$betweenss/kmeans_tissue_list[[k]]$totss)
#   tot_withinss <- c(tot_withinss,  kmeans_tissue_list[[k]]$tot.withinss)
#   w_ss <- c(w_ss,kmeans_tissue_list[[k]]$withinss)
# }
#
# plot(2:31, tot_withinss, type = "l")
# abline(v = 3, lty = 3)
#
#
# #
# tissue_kmeans <- kmeans(tissue_df, centers = 3,  nstart = 10, iter.max = 100)

################################################################################
# Tissue - Hierarchical Clustering
################################################################################
h_tissue <- hclust(dist(tissue_df), method = "ward.D")
plot(h_tissue, main = "Tissue dendrogram",
     xlab = "Observation", ylab = "Height")

# Combine Tissue and Cell dendrograms (for inclusion in thesis)
par(mfrow = c(1,2))
h_tissue %>% as.dendrogram() %>% set("branches_k_color", k = 4) %>%
  set("labels", NA) %>%
  plot(ylab = "Height", main = "Tissue dendrogram")

h %>% as.dendrogram() %>% set("branches_k_color", k = 4) %>%
  set("labels", NA) %>%
  plot(ylab = "Height", main = "Cell dendrogram")


ch_tissue <- cutree(h_tissue, k = 4)
table(ch_tissue)
library(ggrepel)
ggplot(data.frame(x = seq(1:dim(tissue_df)[1]), clust = ch_tissue), aes(x = x, y = clust, col = factor(clust))) + geom_point() +
  geom_text_repel(label = ifelse(ch_tissue == 1, rownames(tissue_df), ""))


h_cluster_tissue_df <- left_join(data.frame(RNA_number_id = names(ch_tissue),
                                            cluster = ch_tissue), bio_source_all)
create_ggplot_heatmap(h_cluster_tissue_df$RNA_number_id, h_cluster_tissue_df$cluster)

################################################################################
# Tissue - UMAP Reduction
################################################################################
h_plots <- create_umap_plots(tissue_df, ch_tissue, 1046, 3)
#k_plots <- create_umap_plots(tissue_df, tissue_kmeans$cluster, 5454, 3)

plots <- h_plots[[1]]; plots2 <- h_plots[[2]]
# Compare plots
split <- length(plots)/2

pdf(here("DR Files", "Tissue UMAP plots h_1.pdf"), 20, 10)
patchwork::wrap_plots(plots[1:split], nrow = 4, ncol = 4)
patchwork::wrap_plots(plots[(split + 1):length(plots)], nrow = 4, ncol = 4)
dev.off()

pdf(here("DR Files", "Tissue UMAP plots h_2.pdf"), 20, 10)
patchwork::wrap_plots(plots2[1:split], nrow = 4, ncol = 4)
patchwork::wrap_plots(plots2[(split + 1):length(plots2)], nrow = 4, ncol = 4)
dev.off()

tissue_umap <- create_umap_single(tissue_df, 35, 0.001, 3, 1046)
ggpairs(data.frame(tissue_umap$layout, "clust" = ch_tissue), columns = 1:3, aes(color = factor(clust)),
        upper = list(continuous = "blank"), diag = list(continuous = "blankDiag"))

saveRDS(tissue_umap$layout, here("RNA Splicing Data", "Tissue UMAP.RDS"))
