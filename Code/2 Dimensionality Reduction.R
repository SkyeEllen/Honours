################################################################################
# Setup
################################################################################
# Libraries
here::i_am("Code/2 Dimensionality Reduction.R")
library(here)         # Locate files relative to this one
library(umap)         # Calculate UMAP embeddings
library(ggplot2)      # Visualisation of UMAP embeddings / heatmaps
library(progress)     # Progress bar for creating umap plots
library(patchwork)    # Plot umap options to pdf in a grid
library(GGally)       # Pairwise plots
library(tidyverse)    # Data manipulation stuff
library(dendextend)   # Changing dendrogram aesthetics
library(RColorBrewer) # Colors

# My additional functions for plotting
source(here("Code", "DR and EA Functions.R"))

# Read in Data
cell_df <- readRDS(here("RNA Splicing Data", "Cell Data.RDS"))
tissue_df <- readRDS(here("RNA Splicing Data", "Tissue Data.RDS"))
bio_source_all <- readRDS(here("RNA Splicing Data", "Bio_source_all.RDS"))
tissue_source <- readRDS(here("RNA Splicing Data", "Tissue source.RDS"))

################################################################################
# Cell - Hierarchical Clustering
################################################################################
d <- dist(cell_df)
h <- hclust(d, method = "ward.D")  # 4 clusters under complete, 2~4 under ward.D
plot(h)
saveRDS(h, here("RNA Splicing Data", "Cell Hierarchical.RDS"))

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

library(rgl)
plot3d(x = u$layout[,1], y = u$layout[,2], z = u$layout[,3], col = ch)

cell_umap <- u$layout
saveRDS(cell_umap, here("RNA Splicing Data", "Cell UMAP.RDS"))

################################################################################
# Tissue - Hierarchical Clustering
################################################################################
h_tissue <- hclust(dist(tissue_df), method = "ward.D")
plot(h_tissue, main = "Tissue dendrogram",
     xlab = "Observation", ylab = "Height")
saveRDS(h_tissue, here("RNA Splicing Data", "Tissue Hierarchical.RDS"))

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

# Check separation in 3D plot
library(rgl)
plot3d(x = tissue_umap$layout[,1], y = tissue_umap$layout[,2], z = tissue_umap$layout[,3], col = ch_tissue)

# Save results
saveRDS(tissue_umap$layout, here("RNA Splicing Data", "Tissue UMAP.RDS"))
