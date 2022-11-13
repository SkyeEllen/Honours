###### Results #######
simulation1 <- readRDS(here("Simulation New Results", "SD1.RDS"))
simulation2 <- readRDS(here("Simulation New Results", "SD2.RDS"))
simulation3 <- readRDS(here("Simulation New Results", "SD3.RDS"))
plot_actual(simulation1)
plot_actual(simulation2)
plot_actual(simulation3)

# Setup stuff
sim_opts <- c("sim1", "sim2", "sim3")
sim_names <- c("Simulation 1", "Simulation 2", "Simulation 3")
chain_opts <- c("", "_3D", "_3D_v2", "_long")
chain_names <-  c("", "3D", "3D_v2", "(Long)")
max_res <- 5
titles <- apply(expand.grid(a = sim_opts, b = chain_opts), 1, FUN = "paste", collapse = "")
names(titles) <- apply(expand.grid(a = sim_names, b = chain_names), 1, FUN = "paste", collapse = " ")

# Cutoff for clusters
epsilon = 0.1

# Pick which simulation
sim_name <- "sim1"; chain_name <- ""; chain_ext = "_"
sim <- simulation1; df <- sim$df

# Get results
res <- list();
for(i in 1:5)
  res[[i]] <- readRDS(here("Simulation New Results", paste(sim_name, "_res", i, chain_name, ".RDS", sep = "")))

# Adjust based on long or not
burn_in = 2500; mcmc_iter = 2000; jumps = 5;

# Convergence test
end = burn_in + mcmc_iter*jumps; start = burn_in + 1; idx = seq(start, end, 5)
mcmc_list <- list();
for(i in 1:max_res) {mcmc_list[[i]] <- mcmc(res[[i]]$total_l_p[idx])}
gelman.diag(mcmc_list)

# Plot of chains
col_list <- RColorBrewer::brewer.pal(5, "Set2")
plot(res[[1]]$total_l_p, col = col_list[1], type = "l", ylim = c(min, 0),
     ylab = expression("log p(x | s, y, "*theta*")"),
     xlab = "Iteration", main = "Simulated Dataset 3")
for(i in 2:max_res) lines(res[[i]]$total_l_p, col = col_list[i])
legend("bottomleft", legend = paste("Chain", 1:5), lty = 1, col = cols)

# Get median K
total_K_list <- list(); for(i in 1:max_res) total_K_list[[i]] <- res[[i]]$total_K[start:end]
lapply(total_K_list, median)

# Convergence on K (not used)
mcmc_list <- list();
for(i in 1:max_res) {mcmc_list[[i]] <- mcmc(res[[i]]$total_K[idx])}
gelman.diag(mcmc_list)

# Plot of K
plot(res[[1]]$total_K, col = col_list[1], type = "l", ylim = c(0, 15))
for(i in 2:5) lines(res[[i]]$total_K, col = col_list[i])

# Distributional Co-Clustering
sj_c <- rbind(read.csv(here("Simulation New Results", paste("Sj_sim", sim_name, "_res", 1, chain_ext, ".txt", sep = "")),
                       header = F));
for(i in 2:5) {
  sj_c <- rbind(read.csv(here("Simulation New Results", paste("Sj_sim", sim_name, "_res", i, chain_ext, ".txt", sep = "")), header = F),
                sj_c)};

# Calculate Posterior co-clustering probability
dist_cc_mat <- coclusteringMatrix(sj_c) + diag(length(res[[1]]$Sj))
dist_order <- order.hclust(hclust(as.dist(1 - dist_cc_mat)))
heatmap(dist_cc_mat[dist_order,dist_order], scale = "none", Rowv = NA, Colv = NA)

# True cluster
true_clust <- true_clusters(sim, "dist")

# Data frame for plotting
h <- hclust(as.dist(1 - dist_cc_mat))
dist_idx <- order(sim$D)
dist_hm_df <- expand_grid(obs1 = 1:dim(dist_cc_mat)[1], obs2 = 1:dim(dist_cc_mat)[1])
dist_hm_df$cc_prob <- as.numeric(dist_cc_mat)
dist_hm_df$true <- as.numeric(true_clust)

# Create plot of comparison
p <- create_post_cc_heatmaps(dist_hm_df, dist_idx,
                             title = paste(names(titles)[titles == paste(sim_name, chain_name, sep = "")],
                                           "Population level"))

# Save above plot
ggsave(paste(sim_name, "_D_heatmaps.pdf"), p)

######## Full MCMC ############
sim_name <- "sim1"; chain_name <- ""
sim <- simulation2; df <- sim$df
sij_c <- rbind(read.csv(here("Simulation New Results", paste("sj_obs_sim", sim_name, "_res", 1, chain_ext, ".txt", sep = "")),
                        header = F));
rij_c <- rbind(read.csv(here("Simulation New Results", paste("Rij_sim", sim_name, "_res", 1, chain_ext, ".txt", sep = "")),
                        header = F))
for(i in 2:5) {
  sij_c <- rbind(read.csv(here("Simulation New Results",
                               paste("sj_obs_sim", sim_name, "_res", i, ext, ".txt",
                                     sep = "")), header = F),
                 sij_c)
  rij_c <- rbind(read.csv(here("Simulation New Results",
                               paste("Rij_sim", sim_name, "_res", 1, ext, ".txt",
                                     sep = "")), header = F),
                 rij_c)
};
combined <- paste(as.numeric(as.matrix(sij_c)), as.numeric(as.matrix(rij_c)))
combined <- matrix(combined, ncol = length(res[[1]]$sj_obs))

# Calculate Posterior co-clustering probability
obs_cc_mat <- coclusteringMatrix(combined) + diag(length(res[[1]]$sj_obs))
obs_order <- order.hclust(hclust(dist(1 - obs_cc_mat)))
heatmap(obs_cc_mat[obs_order,obs_order], scale = "none", Rowv = NA, Colv = NA)

roc_curve2(sim, obs_cc_mat, "Full MCMC Observational")

# 1 - 1
# Sim 1 - None
# Sim 2 - None
# Sim 3 - None

# 0.8 - 0.8
# Sim 1 - 0.01 - 0.46
# Sim 2 - 0.01 0.02
# Sim 3 - 0.03 - 0.11

h <- hclust(as.dist(1 - obs_cc_mat))
plot(h)
epsilon = 0.007
clusters_h <- cutree(h, h = 1 - epsilon)
table(clusters_h, paste(sim$d, sim$O))
heatmap(table(clusters_h, paste(sim$d, sim$O)), scale = "col")
#
cluster_k <- cutree(h, k = length(unique(paste(sim$d, sim$O))))
table(cluster_k, paste(sim$d, sim$O))
heatmap(table(cluster_k, paste(sim$d, sim$O)), scale = "none")
