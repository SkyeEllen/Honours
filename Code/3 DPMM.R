# Application of the DPMM to both datasets & visualisations/analysis of results
################################################################################
# Cell - Dirichlet Process Mixture Model
################################################################################
library(dirichletprocess)
cell_dpmm <- DirichletProcessMvnormal(u$layout)
cell_dpmm <- Fit(cell_dpmm, 1000)
ggplot(data.frame(u$layout, "clust" = cell_dpmm$clusterLabels), aes(x = X1, y = X2, col = factor(clust))) +
  geom_point()

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
