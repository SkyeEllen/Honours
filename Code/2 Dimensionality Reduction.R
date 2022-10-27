################################################################################
# Setup
################################################################################
# Libraries
here::i_am("Code/2 Dimensionality Reduction.R")
library(here)
library(umap)
library(ggplot2)
library(rgl)
library(progress)
library(ggrepel)
library(cluster)
library(patchwork)


# Read in Data
cell_df <- readRDS(here("RNA Splicing Data", "Cell Data.RDS"))
tissue_df <- readRDS(here("RNA Splicing Data", "Tissue Data.RDS"))
tissue_to_cell_sys_conv <- read.csv(here("RNA Splicing Data", "Tissue to Cell Sys.csv"))

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
  #browser()
  # hm_df <- data.frame(clust = clusters, source = source)
  # u_c <- sort(unique(hm_df$clust))
  # u_s <- sort(unique(hm_df$source))
  #
  #
  # hmmat <- matrix(0, nrow = length(u_s), ncol = length(u_c))
  # rownames(hmmat) <- u_s; colnames(hmmat) <- u_c
  # for(i in 1:length(u_s)){
  #   for(j in 1:length(u_c)){
  #     hmmat[i,j] = dim(hm_df %>% filter(clust == u_c[j] & source == u_s[i]))[1]
  #   }
  # }
  hmmat <- table(source, clusters)
  h <- heatmap(hmmat, col = hcl.colors((max(hmmat) + 1), "Reds 2", rev = TRUE))
  legend(x="topleft", legend=c(1:max(hmmat)),
         fill=hcl.colors(max(hmmat), "Reds 2", rev = TRUE))
  return(h)
}

create_ggplot_heatmap <- function(clusters, source){
  hm_df <- expand.grid(cluster = unique(clusters), source = unique(source))
  count <- as.numeric(table(clusters, source))
  hm_df$count <- ifelse(count != 0, count, NA)
  browser()
  ggplot(hm_df, aes(x=cluster, y=source, fill=count)) +
    geom_tile(color="white", size = 0.25) +
    scale_fill_gradient(low="gold", high="darkorchid", na.value="white") + theme_minimal()
}


################################################################################
# Cell Data - K-means clustering
################################################################################
kmeans_cell_list <- list()
for(k in 1:30){
  cat(k)
  kmeans_cell_list[[k]] <- kmeans(cell_df, centers = (k+1),  nstart = 10, iter.max = 100)
}

bet_p <- c(); w_ss <- c(); tot_withinss <- c()


for(k in 1:30) {
  bet_p <- c(bet_p, kmeans_cell_list[[k]]$betweenss/kmeans_cell_list[[k]]$totss)
  tot_withinss <- c(tot_withinss,  kmeans_cell_list[[k]]$tot.withinss)
  w_ss <- c(w_ss,kmeans_cell_list[[k]]$withinss)
}

plot(2:31, tot_withinss, type = "l")
abline(v = 5)


silhouette_score <- function(kmeans, df){
  ss <- silhouette(kmeans$cluster, dist(df))
  mean(ss[, 3])
}

library(cluster)
ss <- c()
for(k in 1:30){
  cat(k)
  ss <- c(ss, silhouette_score(kmeans_cell_list[[k]], cell_df))
}

plot(2:31, type='b', ss, xlab='Number of clusters', ylab='Average Silhouette Scores', frame=FALSE)

saveRDS(kmeans_cell_list, "kmeans_cell_list.RDS")

#
cell_kmeans <- kmeans(cell_df, centers = 4,  nstart = 10, iter.max = 100)
create_ggplot_heatmap(cell_kmeans$cluster, Cell_source)
################################################################################
# Cell - Hierarchical Clustering
################################################################################
d <- dist(cell_df)
h<- hclust(dist(cell_df))
plot(h)
ch <- cutree(h, k = 4)
ggplot(data.frame(x = seq(1:138), ch), aes(x = x, y = ch, col = factor(ch))) + geom_point() +
  geom_text_repel(label = ifelse(ch == 4, rownames(cell_df), ""))


# cell_df <- cell_df[-96,]
# h2 <- hclust(dist(cell_df))
# plot(h2)
# ch <- cutree(h2, 3)

################################################################################
# Cell - PCA Reduction
################################################################################
pca_cell <- prcomp(cell_df, retx = TRUE)
var_explained <- pca_cell$sdev^2/sum(pca_cell$sdev^2)
plot(var_explained)
cum_var_explained <- cumsum(var_explained)
plot(cum_var_explained, ylim = c(0,1))
points(var_explained)
abline(h = 1/138, lty = 2)
abline(h = 0.9, lty = 2)

cutoff_cell <- max(which(var_explained > 1/138))
h_pca <- hclust(dist(pca_cell$x[,1:cutoff_cell]))
plot(h_pca)
ch_hpca <- cutree(h_pca, h = 50)
# Takes 22 components to reach level of explains > than expected (1/n)
# Takes 105 components to reach 90% of variance explained

################################################################################
# Cell - UMAP Reduction
################################################################################
# Create plots
# Takes ~2 minutes to run
h_plots <- create_umap_plots(cell_df, ch, 5454, 3)
k_plots <- create_umap_plots(cell_df, cell_kmeans$cluster, 5454, 3)

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

plots <- k_plots[[1]]; plots2 <- k_plots[[2]]
pdf(here("DR Files", "Cell UMAP plots k_1.pdf"), 20, 10)
patchwork::wrap_plots(plots[1:split], nrow = 4, ncol = 4)
patchwork::wrap_plots(plots[(split + 1):length(plots)], nrow = 4, ncol = 4)
dev.off()

pdf(here("DR Files", "Cell UMAP plots k_2.pdf"), 20, 10)
patchwork::wrap_plots(plots2[1:split], nrow = 4, ncol = 4)
patchwork::wrap_plots(plots2[(split + 1):length(plots2)], nrow = 4, ncol = 4)
dev.off()

# Get umap with selected configuration
u <- create_umap_single(cell_df, 5, 0.1, 4, 5454)
ggpairs(data.frame(u$layout, "clust" = cell_kmeans$cluster), columns = 1:4, aes(color = factor(clust)),
               upper = list(continuous = "blank"), diag = list(continuous = "blankDiag"))


cell_umap <- u$layout
saveRDS(cell_umap, "Cell UMAP.RDS")

################################################################################
# Tissue - K-means Clustering
################################################################################
kmeans_tissue_list <- list()
for(k in 1:30){
  cat(k)
  kmeans_tissue_list[[k]] <- kmeans(tissue_df, centers = (k+1),  nstart = 10, iter.max = 100)
}

bet_p <- c(); w_ss <- c(); tot_withinss <- c()

for(k in 1:30) {
  bet_p <- c(bet_p, kmeans_tissue_list[[k]]$betweenss/kmeans_tissue_list[[k]]$totss)
  tot_withinss <- c(tot_withinss,  kmeans_tissue_list[[k]]$tot.withinss)
  w_ss <- c(w_ss,kmeans_tissue_list[[k]]$withinss)
}

plot(2:31, tot_withinss, type = "l")
abline(v = 3)

saveRDS(kmeans_tissue_list, "kmeans_tissue_list.RDS")

#
tissue_kmeans <- kmeans(tissue_df, centers = 3,  nstart = 10, iter.max = 100)

################################################################################
# Tissue - Hierarchical Clustering
################################################################################
h_tissue <- hclust(dist(tissue_df))
plot(h_tissue)
ch_tissue <- cutree(h_tissue, k = 4)
ggplot(data.frame(x = seq(1:dim(tissue_df)[1]), clust = ch_tissue), aes(x = x, y = clust, col = factor(clust))) + geom_point() +
  geom_text_repel(label = ifelse(ch_tissue == 1, rownames(tissue_df), ""))

create_heatmap(ch_tissue, tissue_to_cell_sys_conv$Cell.System)

################################################################################
# Tissue - PCA Reduction
################################################################################
pca_tissue <- prcomp(tissue_df, retx = TRUE)
var_explained <- pca_tissue$sdev^2/sum(pca_tissue$sdev^2)
plot(var_explained)
cum_var_explained <- cumsum(var_explained)
plot(cum_var_explained, ylim = c(0,1))
points(var_explained)
abline(h = 1/45, lty = 2)
abline(h = 0.9, lty = 2)

cutoff <- max(which(var_explained > 1/45))
cutoff_max <- min(which(cum_var_explained > 0.9))
h_pca_t <- hclust(dist(pca_tissue$x[,1:cutoff]))
plot(h_pca_t)
clust_hpca_t <- cutree(h_pca_t, k = 5)

################################################################################
# Tissue - UMAP Reduction
################################################################################
h_plots <- create_umap_plots(tissue_df, clust_hpca_t, 5454, 3)
k_plots <- create_umap_plots(tissue_df, tissue_kmeans$cluster, 5454, 3)


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

plots <- k_plots[[1]]; plots2 <- k_plots[[2]]
pdf(here("DR Files", "Tissue UMAP plots k_1.pdf"), 20, 10)
patchwork::wrap_plots(plots[1:split], nrow = 4, ncol = 4)
patchwork::wrap_plots(plots[(split + 1):length(plots)], nrow = 4, ncol = 4)
dev.off()

pdf(here("DR Files", "Tissue UMAP plots k_2.pdf"), 20, 10)
patchwork::wrap_plots(plots2[1:split], nrow = 4, ncol = 4)
patchwork::wrap_plots(plots2[(split + 1):length(plots2)], nrow = 4, ncol = 4)
dev.off()

tissue_umap <- create_umap_single(tissue_df, 20, 0.01, 3, 5454)
ggpairs(data.frame(u$layout, "clust" = clust_hpca_t), columns = 1:3, aes(color = factor(clust)),
        upper = list(continuous = "blank"), diag = list(continuous = "blankDiag"))

saveRDS(tissue_umap$layout, "Tissue UMAP.RDS")

################################################################################
# Cell - Dirichlet Process Mixture Model
################################################################################


################################################################################
# Tissue - Dirichlet Process Mixture Model
################################################################################
