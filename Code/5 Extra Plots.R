####### Cell & Tissue UMAP plots
cell_u <- readRDS(here("RNA Splicing Data", "Cell UMAP.RDS"))
tissue_u <- readRDS(here("RNA Splicing Data", "Tissue UMAP.RDS"))

cell_h <- readRDS(here("RNA Splicing Data", "Cell Hierarchical.RDS"))
tissue_h <- readRDS(here("RNA Splicing Data", "Tissue Hierarchical.RDS"))

library(ggplot2)
cell <- ggplot(data.frame(cell_u, clust = cutree(cell_h, 4)),
               aes(x = X1, y = X2, col = factor(clust))) +
  geom_point() + labs(x = "UMAP Embedding 1",
                      y = "UMAP Embedding 2",
                      title = "Cell UMAP",
                      col = "Cluster")

tissue <- ggplot(data.frame(tissue_u, clust = cutree(tissue_h, 4)),
                 aes(x = X1, y = X2, col = factor(clust))) +
  geom_point() + labs(x = "UMAP Embedding 1",
                      y = "UMAP Embedding 2",
                      title = "Tissue UMAP",
                      col = "Cluster")

cell + tissue
ggsave(here("Images", "UMAPs.pdf"), width = 10, height = 5)


####### Cell & Tissue UMAP Histograms
par(mfrow = c(3,2))

# Select bandwidth using Normal Reference Rule (Math5895)
R_hat <- c()
sigma_hat <- c()
n <- aggregate(energy$Energy, by = list(energy$brand), FUN = length)[,2]
min <- c()
for(i in unique(energy$brand)) {
  R_hat <- c(R_hat, IQR(df_u[,i]))
  sigma_hat <- c(sigma_hat, sd(energy$Energy[energy$brand == i]))
  min <- c(min, min(R_hat[length(R_hat)]/1.34, sigma_hat[length(sigma_hat)]))
}
hnr <- 1.06*min*n^(-1/5)

par(mfrow = c(3,2))
titles_u <- paste("UMAP Embedding", 1:3)
for(i in 1:3){
  for(df_u in list(cell_u, tissue_u)){
    print(df_u[1,i])
    hist(df_u[,i], prob = T, col = "lightgrey",
         main = ifelse(i == 1,
                       ifelse(df_u[1,i] == cell_u[1,1],
                              "Cell",
                              "Tissue"),
                       ""),
         xlab = titles_u[i])
    lines(density(df_u[,i]), col = "red")
  }
}


######## Table with number of observations in each system ########

# Application of the mNDP to both datasets & visualisations/analysis of results
here::i_am("Code/4 Multivariate mNDP.R")
library(here)
library(tidyverse)

cell_umap <- readRDS(here("RNA Splicing Data", "Cell UMAP.RDS"))
tissue_umap <- readRDS(here("RNA Splicing Data", "Tissue UMAP.RDS"))
bio_source_all <- readRDS(here("RNA Splicing Data", "Bio_source_all.RDS"))

cell_source <-  readRDS(here("RNA Splicing Data", "Cell source.RDS"))
tissue_source <-  readRDS(here("RNA Splicing Data", "Tissue source.RDS"))

# Set up df for mNDP
cell_source$Bio_source_num <- as.numeric(as.factor(cell_source$Biological_source))
cell_df <- left_join(cell_source,
                     data.frame(cell_umap, RNA_number_id = rownames(cell_umap)))

tissue_source$Cell_sys_num <- as.numeric(as.factor(tissue_source$Cell_system))
tissue_df <- left_join(tissue_source,
                       data.frame(tissue_umap, RNA_number_id = rownames(tissue_umap)))

latex_str <- ""
for(u in names(table(cell_df$Biological_source))){
  syst <- str_sub(u, end = -13)

  latex_str <- paste(latex_str,
                     syst, sep = "")
  latex_str <- paste(latex_str, length(cell_df[cell_df$Biological_source == u,1]),
                     length(tissue_df[tissue_df$Cell_system == syst,1]),
                     sep = "&")
  latex_str <- paste(latex_str, "\\", sep = "")
}

######## Table with number of observations in each system ########
######## Combine levels
library(tidyr)
tissue_df <- readRDS(here("RNA Splicing Data", "Tissue group df.RDS"))
tissue_comb_clust <- data.frame(RNA_number_id = rownames(tissue_df))
tissue_comb_clust$obs <- tissue_obs_clust_df$posterior_clusters
tissue_comb_clust <- left_join(tissue_comb_clust, tissue_source)

tissue_comb_clust <- left_join(tissue_comb_clust,
                               data.frame(Cell_system = tissue_dist_clust_df$System,
                                          Pop_clust = tissue_dist_clust_df$posterior_clusters))
tissue_comb_clust <- left_join(tissue_comb_clust,
                               data.frame(RNA_number_id = rownames(tissue_df),
                                          X1 = tissue_df$X1,
                                          X2 = tissue_df$X2))

for(i in 1:max(tissue_comb_clust$Pop_clust)){
  tissue_comb_clust$obs[tissue_comb_clust$Pop_clust == i] = exclude.empty(tissue_comb_clust$obs[tissue_comb_clust$Pop_clust == i])
}

ggplot(tissue_comb_clust, aes(x = X1, y = X2, col = factor(obs))) + geom_point() +
  facet_wrap(.~Pop_clust)


clust <- paste(tissue_comb_clust$Pop_clust, tissue_comb_clust$obs)
for(i in unique(tissue_comb_clust$Pop_clust)){
  for(j in unique(tissue_comb_clust$obs)){
    cat(paste(unique(tissue_comb_clust$Cell_system[tissue_comb_clust$Pop_clust == i
                                                   & tissue_comb_clust$obs == j]), collapse = ", "), "\n")
  }
}


library(tidyr)
cell_df <- readRDS(here("RNA Splicing Data", "Cell group df.RDS"))
cell_comb_clust <- data.frame(RNA_number_id = rownames(cell_df))
cell_comb_clust$obs <- cell_obs_clust_df$posterior_clusters
cell_comb_clust <- left_join(cell_comb_clust, cell_source)
cell_comb_clust$Biological_source <- str_sub(cell_comb_clust$Biological_source, end = -13)
cell_comb_clust <- left_join(cell_comb_clust,
                             data.frame(Biological_source = cell_dist_clust_df$System,
                                        Pop_clust = cell_dist_clust_df$posterior_clusters))
cell_comb_clust <- left_join(cell_comb_clust,
                             data.frame(RNA_number_id = rownames(cell_df),
                                        X1 = cell_df$X1,
                                        X2 = cell_df$X2))

for(i in 1:max(cell_comb_clust$Pop_clust)){
  cell_comb_clust$obs[cell_comb_clust$Pop_clust == i] = exclude.empty(cell_comb_clust$obs[cell_comb_clust$Pop_clust == i])
}

ggplot(cell_comb_clust, aes(x = X1, y = X2, col = factor(obs))) + geom_point() +
  facet_wrap(.~Pop_clust) +
  geom_point(data = cell_comb_clust[cell_comb_clust$Biological_source == "Mesenchymal Stem",],
             mapping = aes(x = X1, y = X2, col = factor(obs)), col = "black")


clust <- paste(cell_comb_clust$Pop_clust, cell_comb_clust$obs)
for(i in unique(cell_comb_clust$Pop_clust)){
  for(j in unique(cell_comb_clust$obs)){
    cat(i, " & ", j , "&" , paste(unique(cell_comb_clust$Biological_source[cell_comb_clust$Pop_clust == i
                                                       & cell_comb_clust$obs == j]), collapse = ", "), "\\\\\n")
  }
}


write.csv(tissue_comb_clust, here("New For Collaborator", "Final Tissue Clusters.csv"))
write.csv(cell_comb_clust, here("New For Collaborator", "Final Cell Clusters.csv"))

for(i in 1:max(cell_comb_clust$Pop_clust)){
  sub_groups <- cell_dist_clust_df$System[cell_dist_clust_df$posterior_clusters == i]
  idx <- which(str_sub(cell_obs_clust_df$Biological_source, end = -13) %in% sub_groups)
  if(length(sub_groups) > 1){
    h_sub <- hclust(as.dist(1 - cell_obsCCM[idx, idx]))
    plot(h_sub)
    cell_comb_clust$obs[cell_comb_clust$Pop_clust == i] = cutree(h_sub, h = 1 - 0.1)}
  else{
    cell_comb_clust$obs[cell_comb_clust$Pop_clust == i] = 1}

}

ggplot(cell_comb_clust, aes(x = X1, y = X2, col = factor(obs))) + geom_point() +
  facet_wrap(.~Pop_clust)


for(i in 1:max(tissue_comb_clust$Pop_clust)){
  sub_groups <- tissue_dist_clust_df$System[tissue_dist_clust_df$posterior_clusters == i]
  idx <- which(tissue_obs_clust_df$Cell_system %in% sub_groups)
  h_sub <- hclust(as.dist(1 - tissue_obsCCM[idx, idx]))
  plot(h_sub)
  tissue_comb_clust$obs[tissue_comb_clust$Pop_clust == i] = cutree(h_sub, h = 1 - 0.2)
}

ggplot(tissue_comb_clust, aes(x = X1, y = X2, col = factor(obs))) + geom_point() +
  facet_wrap(.~Pop_clust)


clust <- paste(tissue_comb_clust$Pop_clust, tissue_comb_clust$obs)
for(i in unique(tissue_comb_clust$Pop_clust)){
  for(j in unique(tissue_comb_clust$obs)){
    cat(paste(unique(tissue_comb_clust$Cell_system[tissue_comb_clust$Pop_clust == i
                                                   & tissue_comb_clust$obs == j]), collapse = ", "), "\n")
  }
}


