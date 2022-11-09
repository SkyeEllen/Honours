# Check Simulation Results (Clean Version)
# Simulation of datasets and analysis of results
here::i_am("Code/Sim C - Results.R")
library(here)  # For reading files
library(coda)  # For checking convergence with gelman.diag
library(clusternomics) # For getting posterior co-clustering probabilities
library(dendextend) # for order hclust
source(here("Code", "MCMC result functions.R"))
source(here("Code", "Sim A - Functions.R"))


# Choose which simulation to check
sim_name <- "sim3"   # Options = sim1, sim2, sim3
chain_name <- ""     # Options = _3D, _3D_v2, _v2

# Read objects
sim <- readRDS(here("Simulation Results", paste(sim_name, "_small_dataset.RDS", sep ="")))
df <- sim$df

######### General Setup ############
burn_in = 2000; jumps = 5;
sim_opts <- c("sim1", "sim2", "sim3")
sim_names <- c("Simulation 1", "Simulation 2", "Simulation 3")
chain_opts <- c("", "_3D", "_3D_v2", "_long")
chain_names <-  c("", "3D", "3D_v2", "Long")

titles <- apply(expand.grid(a = sim_opts, b = chain_opts), 1, FUN = "paste", collapse = "")
names(titles) <- apply(expand.grid(a = sim_names, b = chain_names), 1, FUN = "paste", collapse = " ")

######### Plot Results #########
res1 <- readRDS(here("Simulation Results",  paste(sim_name, "_res1", chain_name, ".RDS", sep = "")))
res2 <- readRDS(here("Simulation Results",  paste(sim_name, "_res2", chain_name, ".RDS", sep = "")))
res3 <- readRDS(here("Simulation Results",  paste(sim_name, "_res3", chain_name, ".RDS", sep = "")))

true <- res1
true$Sj <- sim$D; true$sj_obs <- sim$d; true$rij <- sim$O
p0 <- plot_results(sim$df, true)
p1 <- plot_results(sim$df, res1)
p2 <- plot_results(sim$df, res2)
p3 <- plot_results(sim$df, res3)

#gridExtra::grid.arrange(p0, p1, p2, p3)
ggsave(here(paste(sim_name, chain_name, "_results_plots.pdf", sep = "")), gridExtra::grid.arrange(p0, p1, p2, p3))

######## Check convergence ##########
start = burn_in + 1; end = length(res1$total_l_p)

ylim <- c(min(res1$total_l_p, res2$total_l_p, res3$total_l_p),
          max(res1$total_l_p, res2$total_l_p, res3$total_l_p))
plot(res1$total_l_p, type = "l", col = "green",
     ylab = expression("log p(x|s, r, "*theta*")"),
     xlab = "Iteration",
     main = paste(names(titles)[titles == paste(sim_name, chain_name, sep = "")],
                  "Probability Chains"),
     ylim = ylim)
lines(res2$total_l_p, col = "blue")
lines(res3$total_l_p, col = "red")
legend("bottomright", legend = paste("Chain", 1:3), col = c("green", "red", "blue"),
       lty = 1)

mcmc.res1 <- mcmc(data= res1$total_l_p[start:end], thin = jumps)
mcmc.res2 <- mcmc(data= res2$total_l_p[start:end], thin = jumps)
mcmc.res3 <- mcmc(data= res3$total_l_p[start:end], thin = jumps)
mcmc_list <- list(mcmc.res1, mcmc.res2, mcmc.res3)

gelman.diag(mcmc_list, confidence = 0.95, autoburnin = F)

# plot(res1$total_K, type = "l", col = "green")
# lines(res2$total_K, col = "blue")
# lines(res3$total_K, col = "red")
#
# # Check convergence on K
# mcmc.res1 <- mcmc(data= res1$total_K[start:end], thin = jumps)
# mcmc.res2 <- mcmc(data= res2$total_K[start:end], thin = jumps)
# mcmc.res3 <- mcmc(data= res3$total_K[start:end], thin = jumps)
# mcmc_list <- list(mcmc.res1, mcmc.res2, mcmc.res3)
#
# gelman.diag(mcmc_list, confidence = 0.95, transform=F, autoburnin=F,
#             multivariate=FALSE)

######## Distributional Co-clustering Heatmaps ##########
# Read Data
sj1 <- read.csv(here("Simulation Results", paste("Sj_sim_", sim_name, "_res1", chain_name, ".txt", sep = "")), header = F)
sj2 <- read.csv(here("Simulation Results", paste("Sj_sim_", sim_name, "_res2", chain_name, ".txt", sep = "")), header = F)
sj3 <- read.csv(here("Simulation Results", paste("Sj_sim_", sim_name, "_res3", chain_name, ".txt", sep = "")), header = F)
sj_c <- rbind(sj1, sj2, sj3); rm(sj1, sj2, sj3)

# Calculate Posterior co-clustering probability
dist_cc_mat <- coclusteringMatrix(sj_c) + diag(dim(sj_c)[2])
dist_order <- order.hclust(hclust(dist(dist_cc_mat)))
heatmap(dist_cc_mat[dist_order,dist_order], scale = "none", Rowv = NA, Colv = NA)

true_clust <- true_clusters(sim, "dist")

h <- hclust(dist(dist_cc_mat))
dist_idx <- order(sim$D)
dist_hm_df <- expand_grid(obs1 = 1:dim(dist_cc_mat)[1], obs2 = 1:dim(dist_cc_mat)[1])
dist_hm_df$cc_prob <- as.numeric(dist_cc_mat)
dist_hm_df$true <- as.numeric(true_clust)

p <- create_post_cc_heatmaps(dist_hm_df, dist_idx,
                        paste(names(titles)[titles == paste(sim_name, chain_name, sep = "")],
                        "Population level"))

ggsave(here("Simulation Images", paste(sim_name, chain_name, "_D_heatmaps.pdf", sep ="")),
       p, width = 10, h = 5)
######## Distributional Cluster Accuracy ##########
J <- length(sim$D)

# Accuracy check - p > 0.5
(cm <- confusion_matrix(sim, dist_cc_mat, 0.5, "dist"))
(misclass_rate <- (cm[1,2] + cm[2,1])/(J*(J-1)/2))
(sensitivity <-  cm[1,1]/sum(cm[,1]))
(specificity <- cm[2,2]/sum(cm[,2]))
cat(round(misclass_rate,3), "&", round(sensitivity,3), "&", round(specificity,3))

# Accuracy check - 50% quantile
(cm <- confusion_matrix(sim, dist_cc_mat, quantile(dist_cc_mat, 0.5), "dist"))
(misclass_rate <- (cm[1,2] + cm[2,1])/(J*(J-1)/2))
(sensitivity <-  cm[1,1]/sum(cm[,1]))
(specificity <- cm[2,2]/sum(cm[,2]))
cat(round(misclass_rate,3), "&", round(sensitivity,3), "&", round(specificity,3))

#
# # ROC curve
# sens_list <- c(); spec_list <- c()
# for(q in seq(0, 1, 0.05)){
#   cm <- confusion_matrix(sim, dist_cc_mat, quantile(dist_cc_mat, q), "dist")
#   sens_list <- c(sens_list, cm[1,1]/sum(cm[,1]))
#   spec_list <- c(spec_list, cm[2,2]/sum(cm[,2]))
# }
#
# plot(1 - spec_list, sens_list, type = "l")


######## Observational Co-clustering Heatmaps ##########
sij1 <- read.csv(here("Simulation Results", paste("sj_obs_sim_", sim_name, "_res1", chain_name, ".txt", sep = "")), header = F)
sij2 <- read.csv(here("Simulation Results", paste("sj_obs_sim_", sim_name, "_res2", chain_name, ".txt", sep = "")), header = F)
sij3 <- read.csv(here("Simulation Results", paste("sj_obs_sim_", sim_name, "_res3", chain_name, ".txt", sep = "")), header = F)
sij_c <- as.numeric(as.matrix(rbind(sij1, sij2, sij3))); rm(sij1, sij2, sij3)

rij1 <- read.csv(here("Simulation Results", paste("Rij_sim_", sim_name, "_res1", chain_name, ".txt", sep = "")), header = F)
rij2 <- read.csv(here("Simulation Results", paste("Rij_sim_", sim_name, "_res2", chain_name, ".txt", sep = "")), header = F)
rij3 <- read.csv(here("Simulation Results", paste("Rij_sim_", sim_name, "_res3", chain_name, ".txt", sep = "")), header = F)
rij_c <- as.numeric(as.matrix(rbind(rij1, rij2, rij3))); rm(rij1, rij2, rij3)

n <- length(sim$d)

combined <- paste(sij_c, rij_c, sep = "")
combined <- matrix(as.numeric(combined), ncol = n, nrow = length(sij_c)/(3*n))
obs_cc_mat <- coclusteringMatrix(combined) + diag(n)

obs_idx <- order.hclust(hclust(dist(true_clusters(sim))))
obs_hm_df <- expand_grid(obs1 = 1:n, obs2 = 1:n)
obs_hm_df$cc_prob <- as.numeric(obs_cc_mat)
obs_hm_df$true <- as.numeric(true_clusters(sim))

p <- create_post_cc_heatmaps(obs_hm_df, obs_idx,
                        paste(names(titles)[titles == paste(sim_name, chain_name, sep = "")],
                              "Observational Level"))
ggsave(here("Simulation Images", paste(sim_name, chain_name, "_O_heatmaps.pdf", sep ="")),
       p, width = 10, h = 5)

######## Observational Cluster Accuracy ##########
(cm <- confusion_matrix(sim, obs_cc_mat, 0.5, "obs"))
(misclass_rate <- (cm[1,2] + cm[2,1])/(n*(n-1)/2))
(sensitivity <-  cm[1,1]/sum(cm[,1]))
(specificity <- cm[2,2]/sum(cm[,2]))
cat(round(misclass_rate,3), "&", round(sensitivity,3), "&", round(specificity,3))

(cm <- confusion_matrix(sim, obs_cc_mat, quantile(obs_cc_mat, 0.5), "obs"))
(misclass_rate <- (cm[1,2] + cm[2,1])/(n*(n-1)/2))
(sensitivity <-  cm[1,1]/sum(cm[,1]))
(specificity <- cm[2,2]/sum(cm[,2]))
cat(round(misclass_rate,3), "&", round(sensitivity,3), "&", round(specificity,3))
#
# sens_list <- c(); spec_list <- c()
# for(q in seq(0, 1, 0.05)){
#   cm <- confusion_matrix(sim, obs_cc_mat, quantile(obs_cc_mat, q), "obs")
#   sens_list <- c(sens_list, cm[1,1]/sum(cm[,1]))
#   spec_list <- c(spec_list, cm[2,2]/sum(cm[,2]))
# }
#
# plot(1 - spec_list, sens_list, type = "l")
obs_idx <- order.hclust(hclust(dist(obs_cc_mat), method = "ward.D"))

plot(hclust(dist(obs_cc_mat), method = "complete"))
ch <- cutree(hclust(dist(obs_cc_mat), method = "complete"), k = 8)
tr <- as.numeric(factor(paste(sim$d, sim$O)))
caret::confusionMatrix(factor(tr), factor(ch))
