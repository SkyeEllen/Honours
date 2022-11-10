################################################################################
# Cell - Dirichlet Process Mixture Model
################################################################################
library(dirichletprocess)
cell_umap <- readRDS(here("RNA Splicing Data", "Cell UMAP.RDS"))
cell_dpmm <- DirichletProcessMvnormal(cell_umap)
cell_dpmm1 <- Fit(cell_dpmm, 1000)
cell_dpmm2 <- Fit(cell_dpmm, 1000)
cell_dpmm3 <- Fit(cell_dpmm, 1000)

saveRDS(cell_dpmm1, here("DPMM Results", "Cell DPMM 1.RDS"))
saveRDS(cell_dpmm2, here("DPMM Results", "Cell DPMM 2.RDS"))
saveRDS(cell_dpmm3, here("DPMM Results", "Cell DPMM 3.RDS"))


ggplot(data.frame(cell_umap, "clust" = cell_dpmm3$clusterLabels), aes(x = X1, y = X2, col = factor(clust))) +
  geom_point()

cell_dpmm_plot_df <- data.frame("Chain1" = cell_dpmm1$clusterLabels,
                                row.names = rownames(cell_umap),
                                "Chain2" = cell_dpmm2$clusterLabels,
                                "Chain3" = cell_dpmm3$clusterLabels)


plot(cell_dpmm1$n)

bio_source <- read.csv(here("RNA Splicing Data", "RNASeq_Atlas_samples.csv"))
match <- paste0("X", bio_source$RNA_number)
Cell_source <- bio_source$Biological_source[match %in% rownames(cell_df)]

create_heatmap(cell_dpmm$clusterLabels, Cell_source)
p <- create_ggplot_heatmap(cell_dpmm$clusterLabels, Cell_source)


################################################################################
# Tissue - Dirichlet Process Mixture Model
################################################################################
tissue_dpmm <- DirichletProcessMvnormal(tissue_umap$layout)
tissue_dpmm1 <- Fit(tissue_dpmm, 1000)
tissue_dpmm2 <- Fit(tissue_dpmm, 1000)
tissue_dpmm3 <- Fit(tissue_dpmm, 1000)


saveRDS(tissue_dpmm1, here("DPMM Results", "Tissue DPMM 1.RDS"))
saveRDS(tissue_dpmm2, here("DPMM Results", "Tissue DPMM 2.RDS"))
saveRDS(tissue_dpmm3, here("DPMM Results", "Tissue DPMM 3.RDS"))


par(mfrow = c(1,3))
plot(tissue_umap$layout[,1], tissue_umap$layout[,2], col = tissue_dpmm1$clusterLabels)
plot(tissue_umap$layout[,1], tissue_umap$layout[,2], col = tissue_dpmm2$clusterLabels)
plot(tissue_umap$layout[,1], tissue_umap$layout[,2], col = tissue_dpmm3$clusterLabels)
par(mfrow = c(1,1))

library(coda)
for(chain in list(tissue_dpmm1$likelihoodChain, tissue_dpmm2$likelihoodChain, tissue_dpmm3$likelihoodChain))
  chain[is.infinite(chain)] = median(chain)

tissue_dpmm1$likelihoodChain[is.infinite(tissue_dpmm1$likelihoodChain)] = median(tissue_dpmm1$likelihoodChain)
tissue_dpmm2$likelihoodChain[is.infinite(tissue_dpmm2$likelihoodChain)] = median(tissue_dpmm2$likelihoodChain)
tissue_dpmm3$likelihoodChain[is.infinite(tissue_dpmm3$likelihoodChain)] = median(tissue_dpmm3$likelihoodChain)

mcmc.res1 <- mcmc(data= tissue_dpmm1$likelihoodChain, thin = 1)
mcmc.res2 <- mcmc(data= tissue_dpmm2$likelihoodChain, thin = 1)
mcmc.res3 <- mcmc(data= tissue_dpmm3$likelihoodChain, thin = 1)
mcmc.list <- list(mcmc.res1, mcmc.res2, mcmc.res3)
gelman.diag(mcmc.list)

ggplot(data.frame(tissue_umap$layout, "clust" = tissue_dpmm1$clusterLabels), aes(x = X1, y = X2, col = factor(clust))) +
  geom_point()

library(clusternomics)
mat <- matrix(unlist(tissue_dpmm1$labelsChain), ncol = length(tissue_dpmm1$labelsChain[[1]]))
mat <- rbind(mat,
             matrix(unlist(tissue_dpmm2$labelsChain), ncol = length(tissue_dpmm2$labelsChain[[1]])),
             matrix(unlist(tissue_dpmm3$labelsChain), ncol = length(tissue_dpmm3$labelsChain[[1]])))
ccm <- coclusteringMatrix(mat) + diag(45)
heatmap(ccm)

clust <- cutree(hclust(dist(ccm), method = "ward.D"), k = 2)

plot(tissue_umap$layout[,1], tissue_umap$layout[,3], col = factor(clust))

create_heatmap(tissue_dpmm$clusterLabels, Tissue_cell_sys)
p <- create_ggplot_heatmap(tissue_dpmm$clusterLabels, Tissue_cell_sys)

################################################################################
# Combined
################################################################################


################################################################################
# Cell - PCA Reduction
################################################################################
# pca_cell <- prcomp(cell_df, retx = TRUE)
# saveRDS(pca_cell, here("RNA Splicing Data", "Cell PCA Full.RDS"))
pca_cell <- readRDS(here("RNA Splicing Data", "Cell PCA Full.RDS"))
var_explained <- pca_cell$sdev^2/sum(pca_cell$sdev^2)
plot(var_explained)
cum_var_explained <- cumsum(var_explained)
plot(cum_var_explained, ylim = c(0,1), type = "l")
lines(var_explained)
abline(h = 1/138, lty = 2)
abline(h = 0.9, lty = 2)
abline(h = 0.05, lty = 2)

plot(cum_var_explained[0:20], ylim = c(0,1), type = "l")
lines(var_explained[0:20])
abline(h = 1/138, lty = 2)
abline(h = 0.9, lty = 2)

cutoff_cell_max <- min(which(cum_var_explained > 0.9))
abline(v = cutoff_cell_max, col = "blue")
cutoff_cell <- max(which(var_explained > 1/138))
#cutoff_cell <- max(which(var_explained > 0.025))

abline(v = cutoff_cell, col = "red")
abline(v = 6, col = "red")

h_pca <- hclust(dist(pca_cell$x[,1:7]))
plot(h_pca)
ch_hpca <- cutree(h_pca, k = 5)
table(ch_hpca)
# Takes 22 components to reach level of explains > than expected (1/n)
# Takes 105 components to reach 90% of variance explained

library(rgl)
plot3d(x = pca_cell$x[,1], y = pca_cell$x[,2], z = pca_cell$x[,4], col = ch)
ggpairs(as.data.frame(pca_cell$x), columns = 1:4, aes(col = factor(ch)))

cell_pca <- pca_cell$x[,1:cutoff_cell]
saveRDS(cell_pca, here("RNA Splicing Data", "Cell PCA.RDS"))

################################################################################
# Cell - DPMM on PCA Reduction
################################################################################

## PCA
cell_dpmm_pca <- DirichletProcessMvnormal(pca_cell$x[,1:3], numInitialClusters = 1)
cell_dpmm_pca1 <- Fit(cell_dpmm_pca, 1000)
cell_dpmm_pca2 <- Fit(cell_dpmm_pca, 1000)
cell_dpmm_pca3 <- Fit(cell_dpmm_pca, 1000)


saveRDS(cell_dpmm_pca1, here("DPMM Results", "PCA Cell DPMM 1.RDS"))
saveRDS(cell_dpmm_pca2, here("DPMM Results", "PCA Cell DPMM 1.RDS"))
saveRDS(cell_dpmm_pca3, here("DPMM Results", "PCA Cell DPMM 1.RDS"))

ggplot(data.frame(pca_cell$x, "clust" = cell_dpmm_pca3$clusterLabels), aes(x = PC1, y = PC2, col = factor(clust))) +
  geom_point()

hm_df <- data.frame()
create_heatmap(cell_dpmm_pca$clusterLabels, Cell_source)
heatmap_cell <- create_ggplot_heatmap(cell_dpmm_pca$clusterLabels, Cell_source)

################################################################################
# Tissue - PCA Reduction
################################################################################
# pca_tissue <- prcomp(tissue_df, retx = TRUE, scale = F)
# saveRDS(pca_tissue, here("RNA Splicing Data", "Tissue PCA Full.RDS"))
pca_tissue <- readRDS(here("RNA Splicing Data", "Tissue PCA Full.RDS"))
var_explained <- pca_tissue$sdev^2/sum(pca_tissue$sdev^2)
plot(var_explained)
cum_var_explained <- cumsum(var_explained)
plot(cum_var_explained, ylim = c(0,1))
points(var_explained)
abline(h = 1/45, lty = 2)
abline(h = 0.9, lty = 2)

cutoff <- max(which(var_explained > 1/45))
abline(v=cutoff)
cutoff_max <- min(which(cum_var_explained > 0.9))
h_pca_t <- hclust(dist(pca_tissue$x[,1:cutoff]), method = "ward.D")
plot(h_pca_t)
clust_hpca_t <- cutree(h_pca_t, k = 6)


# h_clust_tissue_df_pca <- left_join(data.frame(RNA_number_id = names(clust_hpca_t),
#                                               cluster = clust_hpca_t), tissue_source)
# create_ggplot_heatmap(h_clust_tissue_df_pca$Cell_system, h_clust_tissue_df_pca$cluster)

library(rgl)
plot3d(x = pca_tissue$x[,1], y = pca_tissue$x[,2], z = pca_tissue$x[,3], col = clust_hpca_t)

saveRDS(pca_tissue$x[,1:cutoff], here("RNA Splicing Data", "Tissue PCA.RDS"))

################################################################################
# Tissue - DPMM on PCA Reduction
################################################################################

# PCA
tissue_dpmm_pca <- DirichletProcessMvnormal(pca_tissue$x[,1:cutoff])
tissue_dpmm_pca1 <- Fit(tissue_dpmm, 1000)
tissue_dpmm_pca2 <- Fit(tissue_dpmm, 1000)
tissue_dpmm_pca3 <- Fit(tissue_dpmm, 1000)

saveRDS(tissue_dpmm_pca1, here("DPMM Results", "PCA Tissue DPMM 1.RDS"))
saveRDS(tissue_dpmm_pca2, here("DPMM Results", "PCA Tissue DPMM 2.RDS"))
saveRDS(tissue_dpmm_pca3, here("DPMM Results", "PCA Tissue DPMM 3.RDS"))


ggplot(data.frame(pca_tissue$x, "clust" = tissue_dpmm_pca3$clusterLabels), aes(x = PC1, y = PC2, col = factor(clust))) +
  geom_point()
