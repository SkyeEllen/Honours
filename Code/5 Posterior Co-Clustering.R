here::i_am("Code/5 Posterior Co-Clustering.R")
library(here)
library(clusternomics)
library(stringr)
library(tidyverse)
library(dendextend)
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

sj_c <- read.csv(here("mNDP Results", "Sj_sim_cell_res1.txt"), header = F)
for(i in 2:max_res){sj_c <- rbind(sj_c, read.csv(here("mNDP New Results", paste("Sj_sim_cell_res", i, ".txt", sep = "")),  header = F))}

# Posterior co-clustering probabilities
distributionalCCM <-  coclusteringMatrix(sj_c) + diag(23)
colnames(distributionalCCM) <- rownames(distributionalCCM) <- groups
write.csv(distributionalCCM, file = here("New For Collaborator", "Cell Groups - Posterior Co-clustering probability matrix.csv"))

h <- hclust(as.dist(1 - distributionalCCM))
plot(h)
order_idx <- order.hclust(h)
heatmap(distributionalCCM[order_idx, order_idx], scale = "none")

dist_idx <- groups[order.hclust(h)]
dist_hm_df <- expand_grid(obs1 = groups, obs2 = groups)
dist_hm_df$cc_prob <- as.numeric(distributionalCCM)
dist_hm_df$col <- rep(cutree(h, h = 1 - epsilon_dist), length(groups))

# Plot distributional level groups
(p1 <- ggplot(as.data.frame(dist_hm_df), aes(x = factor(obs1, levels =  dist_idx, ordered = T),
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

cell_dist_clust_df <- data.frame(RNA_number_id = rownames(df))
cell_dist_clust_df$chain1 <- res[[1]]$Sj
cell_dist_clust_df$chain2 <- res[[2]]$Sj
cell_dist_clust_df$chain3 <- res[[3]]$Sj
cell_dist_clust_df$posterior_clusters <- cutree(h, h = 1 - epsilon_dist)

write.csv(cell_obs_clust_df, here("New For Collaborator" ,"cell distribution clusters.csv"))

# Observational level
sij_c <- rbind(res[[1]]$sj_obs)
for(i in 2:max_res) sij_c <- rbind(sij_c, res[[i]]$sj_obs)
sij_c <- as.numeric(as.matrix(sij_c));

# Observational level
rij_c <- rbind(res[[1]]$rij)
for(i in 2:max_res) rij_c <- rbind(rij_c, res[[i]]$rij)
rij_c <- as.numeric(as.matrix(rij_c));

# Combine distributional and observational level observations
n <- dim(df)[1]
combined <- paste(sij_c, rij_c, sep = "")
combined <- matrix(as.numeric(combined), ncol = n)

# Calculate posterior co-clustering probabilities
observationalCCM <-  coclusteringMatrix(combined) + diag(n)
colnames(observationalCCM) <- rownames(observationalCCM) <- rownames(df)
write.csv(observationalCCM, file = here("New For Collaborator", "Cell Observations - Posterior Co-clustering probability matrix.csv"))

h <- hclust(as.dist(1 - observationalCCM))
plot(h)
abline(h = 1 - epsilon_obs, col = "red")
obs_idx <- rownames(df)[order.hclust(h)]
obs_hm_df <- expand_grid(RNA_number_id = rownames(df), RNA_number_id2 = rownames(df))
obs_hm_df$col <- rep(cutree(h, h = 1 - epsilon_obs), length(rownames(df)))
obs_hm_df$cc_prob <- as.numeric(observationalCCM)
obs_hm_df <- left_join(obs_hm_df, cell_source)

# Plot of observational posterior co-clustering probabilities
(p2 <- ggplot(as.data.frame(obs_hm_df), aes(x = factor(RNA_number_id, levels =  obs_idx, ordered = T),
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
plots <- p1 + p2
ggsave(here("New For Collaborator", "Cell Heatmaps.pdf"), plots, width = 10, height = 5)

# Data frame with the different clustering results - observational
cell_obs_clust_df <- data.frame(RNA_number_id = rownames(df))
cell_obs_clust_df$chain1 <- paste(res[[1]]$sj_obs, res[[1]]$rij, sep = "_")
cell_obs_clust_df$chain2 <- paste(res[[2]]$sj_obs, res[[2]]$rij, sep = "_")
cell_obs_clust_df$chain3 <- paste(res[[3]]$sj_obs, res[[3]]$rij, sep = "_")
cell_obs_clust_df$posterior_clusters <- cutree(h, h = 1 - 0.01)
write.csv(cell_obs_clust_df, here("New For Collaborator" ,"cell observation clusters.csv"))

cell_obs_clust_df <- left_join(cell_obs_clust_df, cell_source)
heatmap(table(cell_obs_clust_df$Biological_source,
              cell_obs_clust_df$posterior_clusters), scale = "row")

####### Tissue ########
df <- readRDS(here("RNA Splicing Data", "Tissue group df.RDS"))
max_res <- 3
res <- list()
for(i in 1:max_res) res[[i]] <- readRDS(here("mNDP Results", paste("tissue_res", i, ".RDS", sep = "")))

# Read group data
sj_1 <- read.csv(here("mNDP Results", "Sj_sim_tissue_res1.txt"), header = F)
sj_2 <- read.csv(here("mNDP Results", "Sj_sim_tissue_res2.txt"), header = F)
sj_3 <- read.csv(here("mNDP Results", "Sj_sim_tissue_res3.txt"), header = F)
sj_c <- rbind(sj_1, sj_2, sj_3)

# Posterior co-clustering probabilities
tissue_d_ccM <-  coclusteringMatrix(sj_c) + diag(9)
colnames(tissue_d_ccM) <- rownames(tissue_d_ccM) <- tissue_groups
write.csv(tissue_d_ccM, file = here("For Collaborator", "Tissue Groups - Posterior Co-clustering probability matrix.csv"))

h <- hclust(as.dist(1-tissue_d_ccM))
plot(h)
library(dendextend)
heatmap(tissue_d_ccM, scale = "none")

dist_idx <- tissue_groups[order.hclust(h)]
dist_hm_df <- expand_grid(obs1 = tissue_groups, obs2 = tissue_groups)
dist_hm_df$cc_prob <- as.numeric(tissue_d_ccM)
dist_hm_df$col <- rep(cutree(h, h = 1 - epsilon_dist), 9)

# Plot distributional level groups
p1 <- ggplot(as.data.frame(dist_hm_df), aes(x = factor(obs1, levels =  dist_idx, ordered = T),
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
write.csv(tissue_dist_clust_df, here("For Collaborator" ,"Tissue distribution clusters.csv"))

# Observational level
sij1 <- read.csv(here('mNDP Results', "sj_obs_sim_tissue_res1.txt"), header = F)
sij2 <- read.csv(here('mNDP Results', "sj_obs_sim_tissue_res2.txt"), header = F)
sij3 <- read.csv(here('mNDP Results', "sj_obs_sim_tissue_res3.txt"), header = F)
sij_c <- rbind(sij1, sij2, sij3)

rij1 <- read.csv(here('mNDP Results', "Rij_sim_tissue_res1.txt"), header = F)
rij2 <- read.csv(here('mNDP Results', "Rij_sim_tissue_res2.txt"), header = F)
rij3 <- read.csv(here('mNDP Results', "Rij_sim_tissue_res3.txt"), header = F)
rij_c <- rbind(sij1, sij2, sij3)

combined <- c()
for(i in 1:dim(sij1)[1]){
  combined <- rbind(combined, paste(sij1[i,], rij1[i,], sep = "_"))
}

# Observation level heatmap
library(clusternomics)
observationalCCM <-  coclusteringMatrix(combined) + diag(length(rownames(df)))
colnames(observationalCCM) <- rownames(observationalCCM) <- rownames(df)
write.csv(observationalCCM, file = here("New For Collaborator", "Tissue Observations - Posterior Co-clustering probability matrix.csv"))

h <- hclust(as.dist(1 - observationalCCM))
plot(h)
clusters <- cutree(h, h = 1 - 0.5)
obs_idx <- rownames(df)[order.hclust(h)]
obs_hm_df <- expand_grid(RNA_number_id = rownames(df), RNA_number_id2 = rownames(df))
obs_hm_df$cc_prob <- as.numeric(observationalCCM)
obs_hm_df$col <- rep(cutree(h, h = 1 - 0.5), length(rownames(df)))
obs_hm_df$cc_prob_clust <- as.numeric(ifelse(observationalCCM > quantile(observationalCCM, 0.5), observationalCCM, 0))

# Plot of observational posterior co-clustering probabilities
(p2 <- ggplot(as.data.frame(obs_hm_df), aes(x = factor(RNA_number_id, levels =  obs_idx, ordered = T),
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

plots <- p1 + p2
ggsave(here("New For Collaborator", "Tissue Heatmaps.pdf"), plots, width = 10, height = 5)

# Data frame with the different clustering results - observational
tissue_obs_clust_df <- data.frame(RNA_number_id = rownames(df))
tissue_obs_clust_df$chain1 <- paste(res[[1]]$sj_obs, res[[1]]$rij, sep = "_")
tissue_obs_clust_df$chain2 <- paste(res[[2]]$sj_obs, res[[2]]$rij, sep = "_")
tissue_obs_clust_df$chain3 <- paste(res[[3]]$sj_obs, res[[3]]$rij, sep = "_")
tissue_obs_clust_df$posterior_clusters <- cutree(h, h = 1 - 0.4)
write.csv(cell_obs_clust_df, here("New For Collaborator" ,"Tissue clusters.csv"))

tissue_obs_clust_df <- left_join(tissue_obs_clust_df, tissue_source)
heatmap(table(tissue_obs_clust_df$Cell_system, tissue_obs_clust_df$posterior_clusters), scale = "row")
