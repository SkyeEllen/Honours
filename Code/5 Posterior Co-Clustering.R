here::i_am("Code/5 Posterior Co-Clustering.R")
library(here)
library(clusternomics)
library(tidyverse)
library(dendextend)
source(here("Code", "MCMC result functions.R"))
source(here("Code", "Multivariate mNDP Functions.R"))
source(here("Code", "DR and EA Functions.R"))

cell_source <-  readRDS(here("RNA Splicing Data", "Cell source.RDS"))
groups <- str_sub(sort(unique(cell_source$Biological_source)), end = -13)
tissue_source <- readRDS(here("RNA Splicing Data", "Tissue source.RDS"))
tissue_cell_sys <- readRDS(here("RNA Splicing Data", "Tissue Cell Sys.RDS"))
tissue_groups <- sort(unique(tissue_source$Cell_system))[-c(1,6,12,13)]

# Maximum distance for cutoff
epsilon_dist = 0.08
epsilon_obs = 0.05
####### Cell  ######
df <- readRDS(here("RNA Splicing Data", "Cell group df.RDS"))
max_res <- 5
res <- list()
for(i in 1:max_res) res[[i]] <- readRDS(here("mNDP New Results", paste("cell_res", i, ".RDS", sep = "")))

sj_c <- read.csv(here("mNDP New Results", "Sj_sim_cell_res1.txt"), header = F)
for(i in 2:max_res){sj_c <- rbind(sj_c, read.csv(here("mNDP New Results", paste("Sj_sim_cell_res", i, ".txt", sep = "")),  header = F))}

# Posterior co-clustering probabilities
cell_distCCM <-  coclusteringMatrix(sj_c) + diag(23)
colnames(cell_distCCM) <- rownames(cell_distCCM) <- groups
write.csv(cell_distCCM, file = here("New For Collaborator", "Cell Groups - Posterior Co-clustering probability matrix.csv"))

h <- hclust(as.dist(1 - cell_distCCM))
#par(mfrow = c(1,2))
plot(h, ylab = "", sub = "", main = "Cell System Level Dendrogram", xlab = "")
order_idx <- order.hclust(h)


dist_idx <- groups[order.hclust(h)]
cell_dist_hm_df <- expand_grid(obs1 = groups, obs2 = groups)
cell_dist_hm_df$cc_prob <- as.numeric(cell_distCCM)
cell_dist_hm_df$col <- rep(cutree(h, h = 1 - epsilon_dist), length(groups))

cell_dist_clust_df <- data.frame(System = groups)
cell_dist_clust_df$chain1 <- res[[1]]$Sj
cell_dist_clust_df$chain2 <- res[[2]]$Sj
cell_dist_clust_df$chain3 <- res[[3]]$Sj
cell_dist_clust_df$chain4 <- res[[4]]$Sj
cell_dist_clust_df$chain5 <- res[[5]]$Sj
cell_dist_clust_df$posterior_clusters <- cutree(h, h = 1 - epsilon_dist)
write.csv(cell_dist_clust_df, here("New For Collaborator" ,"cell distribution clusters.csv"))

(cell_hm <- create_ggplot_heatmap(cell_dist_clust_df$posterior_clusters, cell_dist_clust_df$System))

# Observational level distribution allocation
sij_c <- rbind(read.csv(here("mNDP New Results", "sj_obs_sim_cell_res1.txt"), header = F))
for(i in 2:max_res) sij_c <- rbind(sij_c, read.csv(here("mNDP New Results", paste("sj_obs_sim_cell_res", i, ".txt", sep = "")), header = F))
sij_c <- as.numeric(as.matrix(sij_c));

# Observational level observation allocation
rij_c <- rbind(read.csv(here("mNDP New Results", "Rij_sim_cell_res1.txt"), header = F))
for(i in 2:max_res) rij_c <- rbind(rij_c, read.csv(here("mNDP New Results", paste("Rij_sim_cell_res", i, ".txt", sep = "")), header = F))
rij_c <- as.numeric(as.matrix(rij_c));

# Combine distributional and observational level observations
n <- dim(df)[1]
combined <- paste(sij_c, rij_c, sep = "")
combined <- matrix(as.numeric(combined), ncol = n)

# Calculate posterior co-clustering probabilities
cell_obsCCM <-  coclusteringMatrix(combined) + diag(n)
colnames(cell_obsCCM) <- rownames(cell_obsCCM) <- rownames(df)
write.csv(cell_obsCCM, file = here("New For Collaborator", "Cell Observations - Posterior Co-clustering probability matrix.csv"))

h <- hclust(as.dist(1 - cell_obsCCM))
plot(h, ylab = "", sub = "", main = "Cell Observation Level Dendrogram", xlab = "")
abline(h = 1 - epsilon_obs, col = "red")
obs_idx <- rownames(df)[order.hclust(h)]
cell_obs_hm_df <- expand_grid(RNA_number_id = rownames(df), RNA_number_id2 = rownames(df))
cell_obs_hm_df$col <- rep(cutree(h, h = 1 - epsilon_obs), length(rownames(df)))
cell_obs_hm_df$cc_prob <- as.numeric(cell_obsCCM)
cell_obs_hm_df <- left_join(cell_obs_hm_df, cell_source)


# Heatmaps
# Plot distributional level groups
(p1 <- ggplot(as.data.frame(cell_dist_hm_df), aes(x = factor(obs1, levels =  dist_idx, ordered = T),
                                             y = factor(obs2, levels =  dist_idx, ordered = F),
                                             fill = cc_prob,
                                             col = factor(col))) +
    geom_tile(size = 0.3) + theme_bw() +
    theme(axis.ticks = element_blank(), axis.title = element_blank(),
          axis.text.x = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) +
    scale_fill_gradient(low = "#ffffff", high = "#000000") +
    labs(fill = "Co-clustering\n  probability", col = "Cluster", title = "Cell - System Level"))


# Plot of observational posterior co-clustering probabilities
(p2 <- ggplot(as.data.frame(cell_obs_hm_df), aes(x = factor(RNA_number_id, levels =  obs_idx, ordered = T),
                                            y = factor(RNA_number_id2, levels =  obs_idx, ordered = F),
                                            fill = cc_prob, col = factor(col))) +
    geom_tile() + theme_bw() +
    theme(axis.ticks = element_blank(), axis.title = element_blank(),
          axis.text.x = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(), axis.text.y = element_blank()) +
    scale_fill_gradient(low = "#ffffff", high = "#000000") +
    labs(fill = "Co-clustering\n  probability", col = "Cluster", title = "Cell - Observational Level"))
plots <- gridExtra::grid.arrange(p1, p2, ncol = 2)
ggsave(here("New For Collaborator", "Cell Heatmaps.pdf"), plots, width = 10, height = 5)

# Data frame with the different clustering results - observational
cell_obs_clust_df <- data.frame(RNA_number_id = rownames(df))
cell_obs_clust_df$chain1 <- paste(res[[1]]$sj_obs, res[[1]]$rij, sep = "_")
cell_obs_clust_df$chain2 <- paste(res[[2]]$sj_obs, res[[2]]$rij, sep = "_")
cell_obs_clust_df$chain3 <- paste(res[[3]]$sj_obs, res[[3]]$rij, sep = "_")
cell_obs_clust_df$chain4 <- paste(res[[4]]$sj_obs, res[[4]]$rij, sep = "_")
cell_obs_clust_df$chain5 <- paste(res[[5]]$sj_obs, res[[5]]$rij, sep = "_")
cell_obs_clust_df$posterior_clusters <- cutree(h, h = 1 - epsilon_obs)
write.csv(cell_obs_clust_df, here("New For Collaborator" ,"cell observation clusters.csv"))

par(mfrow = c(1,1))
cell_obs_clust_df <- left_join(cell_obs_clust_df, cell_source)
heatmap(table(cell_obs_clust_df$Biological_source,
              cell_obs_clust_df$posterior_clusters), scale = "row")

(cell_hm <- create_ggplot_heatmap(cell_obs_clust_df$posterior_clusters, cell_obs_clust_df$Biological_source))

####### Tissue ########
df <- readRDS(here("RNA Splicing Data", "Tissue group df.RDS"))
max_res <- 5
res <- list()
for(i in 1:max_res) res[[i]] <- readRDS(here("mNDP New Results", paste("tissue_res", i, ".RDS", sep = "")))

# Read group data
sj_c <- read.csv(here("mNDP New Results", "Sj_sim_tissue_res1.txt"), header = F)
for(i in 2:max_res){sj_c <- rbind(sj_c, read.csv(here("mNDP New Results", paste("Sj_sim_tissue_res", i, ".txt", sep = "")),  header = F))}

# Posterior co-clustering probabilities
tissue_d_ccM <-  coclusteringMatrix(sj_c) + diag(9)
colnames(tissue_d_ccM) <- rownames(tissue_d_ccM) <- tissue_groups
write.csv(tissue_d_ccM, file = here("New For Collaborator", "Tissue Groups - Posterior Co-clustering probability matrix.csv"))

par(mfrow = c(1,2))
h <- hclust(as.dist(1-tissue_d_ccM))
plot(h, ylab = "", sub = "", main = "Tissue System Level Dendrogram", xlab = "")
library(dendextend)
#heatmap(tissue_d_ccM, scale = "none")

dist_idx <- tissue_groups[order.hclust(h)]
tissue_dist_hm_df <- expand_grid(obs1 = tissue_groups, obs2 = tissue_groups)
tissue_dist_hm_df$cc_prob <- as.numeric(tissue_d_ccM)
tissue_dist_hm_df$col <- rep(cutree(h, h = 1 - epsilon_dist), 9)

# Plot distributional level groups
p1 <- ggplot(as.data.frame(tissue_dist_hm_df), aes(x = factor(obs1, levels =  dist_idx, ordered = T),
                                            y = factor(obs2, levels =  dist_idx, ordered = F),
                                            fill = cc_prob,
                                            col = factor(col))) +
  geom_tile(size = 0.3) + theme_bw() +
  theme(axis.ticks = element_blank(), axis.title = element_blank(),
        axis.text.x = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  scale_fill_gradient(low = "#ffffff", high = "#000000") +
  labs(fill = "Co-clustering\n  probability", col = "Cluster", title = "System Level")

tissue_dist_clust_df <- data.frame(System = tissue_groups)
tissue_dist_clust_df$posterior_clusters <- cutree(h, h = 1 - epsilon_dist)
tissue_dist_clust_df$chain1 <- res[[1]]$Sj
tissue_dist_clust_df$chain2 <- res[[2]]$Sj
tissue_dist_clust_df$chain3 <- res[[3]]$Sj
tissue_dist_clust_df$chain4 <- res[[4]]$Sj
tissue_dist_clust_df$chain5 <- res[[5]]$Sj
write.csv(tissue_dist_clust_df, here("New For Collaborator" ,"Tissue distribution clusters.csv"))

# Observational level
# Observational level distribution allocation
sij_c <- rbind(read.csv(here("mNDP New Results", "sj_obs_sim_tissue_res1.txt"), header = F))
for(i in 2:max_res) sij_c <- rbind(sij_c, read.csv(here("mNDP New Results", paste("sj_obs_sim_tissue_res", i, ".txt", sep = "")), header = F))
sij_c <- as.numeric(as.matrix(sij_c));

# Observational level observation allocation
rij_c <- rbind(read.csv(here("mNDP New Results", "Rij_sim_tissue_res1.txt"), header = F))
for(i in 2:max_res) rij_c <- rbind(rij_c, read.csv(here("mNDP New Results", paste("Rij_sim_tissue_res", i, ".txt", sep = "")), header = F))
rij_c <- as.numeric(as.matrix(rij_c));

# Combine distributional and observational level observations
n <- dim(df)[1]
combined <- paste(sij_c, rij_c, sep = "")
combined <- matrix(as.numeric(combined), ncol = n)

# Observation level heatmap
library(clusternomics)
tissue_obsCCM <-  coclusteringMatrix(combined) + diag(length(rownames(df)))
colnames(tissue_obsCCM) <- rownames(tissue_obsCCM) <- rownames(df)
write.csv(tissue_obsCCM, file = here("New For Collaborator", "Tissue Observations - Posterior Co-clustering probability matrix.csv"))

h <- hclust(as.dist(1 - tissue_obsCCM))
plot(h, ylab = "", sub = "", main = "Tissue Observation Level Dendrogram", xlab = "")
obs_idx <- rownames(df)[order.hclust(h)]
tissue_obs_hm_df <- expand_grid(RNA_number_id = rownames(df), RNA_number_id2 = rownames(df))
tissue_obs_hm_df$cc_prob <- as.numeric(tissue_obsCCM)
tissue_obs_hm_df$col <- rep(cutree(h, h = 1 - epsilon_obs), length(rownames(df)))

# Plot of observational posterior co-clustering probabilities
(p2 <- ggplot(as.data.frame(tissue_obs_hm_df), aes(x = factor(RNA_number_id, levels =  obs_idx, ordered = T),
                                            y = factor(RNA_number_id2, levels =  obs_idx, ordered = F),
                                            fill = cc_prob,
                                            col = factor(col))) +
    geom_tile(size = 0.3) + theme_bw() +
    theme(axis.ticks = element_blank(), axis.title = element_blank(),
          axis.text.x = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(), axis.text.y = element_blank()) +
    scale_fill_gradient(low = "#ffffff", high = "#000000") +
    labs(fill = "Co-clustering\n  probability", col = "Cluster", title = "Tissue - Observational Level",
         col = "Cluster"))

plots <- gridExtra::grid.arrange(p1, p2, ncol = 2)
ggsave(here("New For Collaborator", "Tissue Heatmaps.pdf"), plots, width = 10, height = 5)

# Data frame with the different clustering results - observational
tissue_obs_clust_df <- data.frame(RNA_number_id = rownames(df))
tissue_obs_clust_df$chain1 <- paste(res[[1]]$sj_obs, res[[1]]$rij, sep = "_")
tissue_obs_clust_df$chain2 <- paste(res[[2]]$sj_obs, res[[2]]$rij, sep = "_")
tissue_obs_clust_df$chain3 <- paste(res[[3]]$sj_obs, res[[3]]$rij, sep = "_")
tissue_obs_clust_df$chain4 <- paste(res[[4]]$sj_obs, res[[4]]$rij, sep = "_")
tissue_obs_clust_df$chain5 <- paste(res[[5]]$sj_obs, res[[5]]$rij, sep = "_")
tissue_obs_clust_df$posterior_clusters <- cutree(h, h = 1 - epsilon_obs)
write.csv(tissue_obs_clust_df, here("New For Collaborator" ,"Tissue clusters.csv"))

tissue_obs_clust_df <- left_join(tissue_obs_clust_df, tissue_source)
heatmap(table(tissue_obs_clust_df$Cell_system, tissue_obs_clust_df$posterior_clusters), scale = "row")

######## Combine levels
library(tidyr)
tissue_df <- readRDS(here("RNA Splicing Data", "Tissue group df.RDS"))
tissue_comb_clust <- data.frame(RNA_number_id = rownames(df))
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
  facet_wrap(.~Pop_clust)

write.csv(tissue_comb_clust, here("New For Collaborator", "Final Tissue Clusters.csv"))
write.csv(cell_comb_clust, here("New For Collaborator", "Final Cell Clusters.csv"))
