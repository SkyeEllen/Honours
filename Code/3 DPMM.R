here::i_am("Code/3 DPMM.R")
library(here) # File location
library(ggplot2) # Data manipulation and plotting
library(dirichletprocess) # Run the DPMM
library(coda)
source(here("Code", "DR and EA Functions.R"))

################################################################################
# Cell - Dirichlet Process Mixture Model
################################################################################
cell_umap <- readRDS(here("RNA Splicing Data", "Cell UMAP.RDS"))
cell_dpmm <- DirichletProcessMvnormal(cell_umap, numInitialClusters = 1)
cell_dpmm1 <- Fit(cell_dpmm, 10000)
cell_dpmm <- DirichletProcessMvnormal(cell_umap, numInitialClusters = 5)
cell_dpmm2 <- Fit(cell_dpmm, 10000)
cell_dpmm <- DirichletProcessMvnormal(cell_umap, numInitialClusters = 10)
cell_dpmm3 <- Fit(cell_dpmm, 10000)

saveRDS(cell_dpmm1, here("DPMM Results", "Cell DPMM 1.RDS"))
saveRDS(cell_dpmm2, here("DPMM Results", "Cell DPMM 2.RDS"))
saveRDS(cell_dpmm3, here("DPMM Results", "Cell DPMM 3.RDS"))


################################################################################
# Cell - DPMM Results
################################################################################
cell_dpmm1 <- readRDS(here("DPMM Results", "Cell DPMM 1.RDS"))
cell_dpmm2 <- readRDS(here("DPMM Results", "Cell DPMM 2.RDS"))
cell_dpmm3 <- readRDS(here("DPMM Results", "Cell DPMM 3.RDS"))
cell_umap <- readRDS(here("RNA Splicing Data", "Cell UMAP.RDS"))

cell_dpmm_plot_df <- data.frame("RNA_number_id" = rownames(cell_umap),
                                "Chain1" = cell_dpmm1$clusterLabels,
                                "Chain2" = cell_dpmm2$clusterLabels,
                                "Chain3" = cell_dpmm3$clusterLabels,
                                cell_umap)

cp1 <- ggplot(data.frame(cell_dpmm_plot_df),
       aes(x = X1, y = X2, col = factor(Chain1))) +
  geom_point()+ labs(x = "UMAP Embedding 1", y = "UMAP Embedding 2",
                     col = "Cluster", title = "Chain 1")

cp2 <- ggplot(data.frame(cell_dpmm_plot_df),
             aes(x = X1, y = X2, col = factor(Chain2))) +
  geom_point() + labs(x = "UMAP Embedding 1", y = "",
                      col = "Cluster", title = "Chain 2")

cp3 <- ggplot(data.frame(cell_dpmm_plot_df),
             aes(x = X1, y = X2, col = factor(Chain3))) +
  geom_point() + labs(x = "UMAP Embedding 1", y = "",
                      col = "Cluster", title = "Chain 3")

cp1 + cp2 + cp3
ggsave(here("DPMM Results", "DPMM Cell clusters.png"), width = 18, height = 6)

# Biological Source Heatmap
cell_source <- readRDS(here("RNA Splicing Data", "Cell source.RDS"))
cell_dpmm_plot_df <- left_join(cell_dpmm_plot_df, cell_source)
p <- create_ggplot_heatmap(cell_dpmm_plot_df$Chain1, cell_dpmm_plot_df$Biological_source)
p

p <- create_ggplot_heatmap(cell_dpmm_plot_df$Chain2, cell_dpmm_plot_df$Biological_source)
p

#
library(coda)

cell_dpmm1$likelihoodChain1 = cell_dpmm1$likelihoodChain
cell_dpmm2$likelihoodChain1 = cell_dpmm2$likelihoodChain
cell_dpmm3$likelihoodChain1 = cell_dpmm3$likelihoodChain
cell_dpmm1$likelihoodChain1[is.infinite(cell_dpmm1$likelihoodChain)] = median(cell_dpmm1$likelihoodChain[!is.infinite(cell_dpmm1$likelihoodChain)])
cell_dpmm2$likelihoodChain1[is.infinite(cell_dpmm2$likelihoodChain)] = median(cell_dpmm2$likelihoodChain[!is.infinite(cell_dpmm2$likelihoodChain)])
cell_dpmm3$likelihoodChain1[is.infinite(cell_dpmm3$likelihoodChain)] = median(cell_dpmm3$likelihoodChain[!is.infinite(cell_dpmm3$likelihoodChain)])

ggplot(data.frame(idx = 1:5000,
                  chain1 = cell_dpmm1$likelihoodChain1,
                  chain2 = cell_dpmm2$likelihoodChain1,
                  chain3 = cell_dpmm3$likelihoodChain1),
                  aes(x = idx)) +
  geom_line(aes(y = chain1), alpha = 0.5, col = "red")+
  geom_line(aes(y = chain2), alpha = 0.5, col = "green")+
  geom_line(aes(y = chain3), alpha = 0.5, col = "blue")

mcmc.res1 <- mcmc(data= cell_dpmm1$likelihoodChain1, thin = 1)
mcmc.res2 <- mcmc(data= cell_dpmm2$likelihoodChain1, thin = 1)
mcmc.res3 <- mcmc(data= cell_dpmm3$likelihoodChain1, thin = 1)
mcmc.list <- list(mcmc.res1, mcmc.res2, mcmc.res3)
gelman.diag(mcmc.list, autoburnin = T, transform = F)

# K
cellK1_list <- unlist(lapply(cell_dpmm1$labelsChain, FUN = "max"))
cellK2_list <- unlist(lapply(cell_dpmm2$labelsChain, FUN = "max"))
cellK3_list <- unlist(lapply(cell_dpmm3$labelsChain, FUN = "max"))
gelman.diag(list(mcmc(cellK1_list), mcmc(cellK2_list), mcmc(cellK3_list)))

idx <- seq(5,length(cellK1_list), by = 5)
plot(cellK1_list[idx], col = "green", type = "l")
lines(cellK2_list[idx], col = "blue")
lines(cellK3_list[idx], col = "red")

# Post cc (given different clusterings at the end)
library(clusternomics)
l = length(unlist(cell_dpmm1$labelsChain))

# Use final clusters
post_cc_cell <- coclusteringMatrix(rbind(cell_dpmm_plot_df$Chain1,
                                         cell_dpmm_plot_df$Chain2,
                                          cell_dpmm_plot_df$Chain3))
heatmap(ifelse(post_cc_cell > 0.5, post_cc_cell, 0),
        scale = "none")
h <- hclust(dist(post_cc_cell))

plot(h)
clust <- cutree(h, k = 3)
create_ggplot_heatmap(clust, cell_dpmm_plot_df$Biological_source)


################################################################################
# Tissue - Dirichlet Process Mixture Model
################################################################################
tissue_umap <- readRDS(here("RNA Splicing Data", "Tissue UMAP.RDS"))
tissue_dpmm <- DirichletProcessMvnormal(tissue_umap, numInitialClusters = 1)
tissue_dpmm1 <- Fit(tissue_dpmm, 5000)
tissue_dpmm <- DirichletProcessMvnormal(tissue_umap, numInitialClusters = 5)
tissue_dpmm2 <- Fit(tissue_dpmm, 5000)
tissue_dpmm <- DirichletProcessMvnormal(tissue_umap, numInitialClusters = 10)
tissue_dpmm3 <- Fit(tissue_dpmm, 5000)


saveRDS(tissue_dpmm1, here("DPMM Results", "Tissue DPMM 1.RDS"))
saveRDS(tissue_dpmm2, here("DPMM Results", "Tissue DPMM 2.RDS"))
saveRDS(tissue_dpmm3, here("DPMM Results", "Tissue DPMM 3.RDS"))

################################################################################
# Tissue - DPMM Results
################################################################################
tissue_umap <- readRDS(here("RNA Splicing Data", "Tissue UMAP.RDS"))
tissue_dpmm1 <- readRDS(here("DPMM Results", "Tissue DPMM 1.RDS"))
tissue_dpmm2 <- readRDS(here("DPMM Results", "Tissue DPMM 2.RDS"))
tissue_dpmm3 <- readRDS(here("DPMM Results", "Tissue DPMM 3.RDS"))

tissue_plot_df <- data.frame("RNA_number_id" = rownames(tissue_umap),
                                "Chain1" = tissue_dpmm1$clusterLabels,
                                "Chain2" = tissue_dpmm2$clusterLabels,
                                "Chain3" = tissue_dpmm3$clusterLabels,
                                tissue_umap)

tp1 <- ggplot(data.frame(tissue_plot_df),
             aes(x = X1, y = X2, col = factor(Chain1))) +
  geom_point()+ labs(x = "UMAP Embedding 1", y = "UMAP Embedding 2",
                     col = "Cluster", title = "Chain 1")

tp2 <- ggplot(data.frame(tissue_plot_df),
             aes(x = X1, y = X2, col = factor(Chain2))) +
  geom_point() + labs(x = "UMAP Embedding 1", y = "",
                      col = "Cluster", title = "Chain 2")

tp3 <- ggplot(data.frame(tissue_plot_df),
             aes(x = X1, y = X2, col = factor(Chain3))) +
  geom_point() + labs(x = "UMAP Embedding 1", y = "",
                      col = "Cluster", title = "Chain 3")

tp1 + tp2 + tp3
ggsave(here("DPMM Results", "DPMM Tissue clusters.png"), width = 18, height = 6)


gridExtra::grid.arrange(cp1, tp1, cp2, tp2, cp3, tp3, ncol = 2)
#cp1 + tp1 + cp2 + tp2 + cp3 + tp3
ggsave(here("DPMM Results", "DPMM All clusters.png"), width = 10, height = 15)


library(coda)

plot(tissue_dpmm1$likelihoodChain, col = "red", type = "l")
lines(tissue_dpmm2$likelihoodChain, col = "green")
lines(tissue_dpmm3$likelihoodChain, col = "blue")

tissue_dpmm1$likelihoodChain[is.infinite(tissue_dpmm1$likelihoodChain)] = min(tissue_dpmm1$likelihoodChain[!is.infinite(tissue_dpmm1$likelihoodChain)])
tissue_dpmm2$likelihoodChain[is.infinite(tissue_dpmm2$likelihoodChain)] = min(tissue_dpmm1$likelihoodChain[!is.infinite(tissue_dpmm1$likelihoodChain)])
tissue_dpmm3$likelihoodChain[is.infinite(tissue_dpmm3$likelihoodChain)] = min(tissue_dpmm1$likelihoodChain[!is.infinite(tissue_dpmm1$likelihoodChain)])

mcmc.res1 <- mcmc(data= tissue_dpmm1$likelihoodChain, thin = 1)
mcmc.res2 <- mcmc(data= tissue_dpmm2$likelihoodChain, thin = 1)
mcmc.res3 <- mcmc(data= tissue_dpmm3$likelihoodChain, thin = 1)
mcmc.list <- list(mcmc.res1, mcmc.res2, mcmc.res3)
gelman.diag(mcmc.list)



# Biological Source Heatmap
tissue_source <- readRDS(here("RNA Splicing Data", "Tissue source.RDS"))
tissue_hm_df <- left_join(tissue_plot_df, tissue_source)
#create_heatmap(tissue_dpmm$clusterLabels, Tissue_cell_sys)
(p <- create_ggplot_heatmap(tissue_hm_df$Chain1, tissue_hm_df$Cell_system))
ggsave(here("Images", "Tissue DPMM heatmap.pdf"))

# K
K1_list <- unlist(lapply(tissue_dpmm1$labelsChain, FUN = "max"))
K2_list <- unlist(lapply(tissue_dpmm2$labelsChain, FUN = "max"))
K3_list <- unlist(lapply(tissue_dpmm3$labelsChain, FUN = "max"))
gelman.diag(list(mcmc(K1_list), mcmc(K2_list), mcmc(K3_list)))


ggplot(data.frame(idx = 1:length(K1_list),
                  chain1 = K1_list,
                  chain2 = K2_list,
                  chain3 = K3_list),
       aes(x = idx, y = chain1)) + geom_line(aes(col = "Chain 1")) +
  geom_line(aes(y = chain2, col = "Chain 2")) +
  geom_line(aes(y = chain3, col = "Chain 3")) + theme_minimal() +
  scale_y_continuous(breaks = seq(2,10, by = 2)) +
  labs(x = "Iteration", y = "K")


plot(K1_list, lty = 1, type = "l")
lines(K2_list, lty = 2)
lines(K3_list, lty = 3)



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
