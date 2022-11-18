here::i_am("Code/3 DPMM.R")
library(here) # File location
library(ggplot2) # Data manipulation and plotting
library(dirichletprocess) # Run the DPMM
library(coda)
library(clusternomics)
library(tidyverse)
library(RColorBrewer)

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

###### Post CC #####
# All converged to same iterations - just use this for labelling purposes
post_cc_cell <- coclusteringMatrix(rbind(cell_dpmm_plot_df$Chain1,
                                         cell_dpmm_plot_df$Chain2,
                                          cell_dpmm_plot_df$Chain3))
h <- hclust(as.dist(1-post_cc_cell))

plot(h)
clust <- cutree(h, h = 1 - 0.08)
cell_hm <- create_ggplot_heatmap(clust, str_sub(cell_dpmm_plot_df$Biological_source, end = -13))
cell_p <- plot_dirichletprocess_multivariate(cell_dpmm3) +
  labs(x = "UMAP Embedding 1", y = "UMAP Embedding 2") + scale_color_brewer(palette = "Dark2")


################################################################################
# Tissue - Dirichlet Process Mixture Model
################################################################################
tissue_umap <- readRDS(here("RNA Splicing Data", "Tissue UMAP.RDS"))
tissue_dpmm <- DirichletProcessMvnormal(tissue_umap, numInitialClusters = 1)
tissue_dpmm1 <- Fit(tissue_dpmm, 10000)
tissue_dpmm <- DirichletProcessMvnormal(tissue_umap, numInitialClusters = 5)
tissue_dpmm2 <- Fit(tissue_dpmm, 10000)
tissue_dpmm <- DirichletProcessMvnormal(tissue_umap, numInitialClusters = 10)
tissue_dpmm3 <- Fit(tissue_dpmm, 10000)


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
ggsave(here("DPMM Results", "DPMM All clusters.png"), width = 10, height = 15)

tissue_p <- plot_dirichletprocess_multivariate(tissue_dpmm1) + #theme_minimal() +
  labs(x = "UMAP Embedding 1", y = "UMAP Embedding 2") +
  scale_color_brewer(palette = "Dark2")
cell_p + tissue_p
ggsave(here("Images", "DPMM Posterior Clusters.pdf"), cell_p + tissue_p, width = 10, height = 5)

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
(tissue_hm <- create_ggplot_heatmap(tissue_hm_df$Chain1, tissue_hm_df$Cell_system))
ggsave(here("Images", "Tissue DPMM heatmap.pdf"))

ggsave(here("Images", "DPMM heatmaps.pdf"), cell_hm + tissue_hm, width = 10, height = 5)

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
