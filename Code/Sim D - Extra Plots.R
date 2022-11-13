here::i_am("Code/Sim D - Extra Plots.R")
library(here)
library(clusternomics)
source(here("Code", "MCMC result functions.R"))

######### General Setup ############
simulation1 <- readRDS(here("Simulation New Results", "SD1.RDS"))
simulation2 <- readRDS(here("Simulation New Results", "SD2.RDS"))
simulation3 <- readRDS(here("Simulation New Results", "SD3.RDS"))

# Choose which simulation to check
# Read objects
burn_in = 2000; jumps = 5;
sim_opts <- c("sim1", "sim2", "sim3")
sim_names <- c("Simulation 1", "Simulation 2", "Simulation 3")
chain_opts <- c("", "_3D", "_3D_v2", "_long")
chain_names <-  c("", "3D", "3D_v2", "(Long)")

titles <- apply(expand.grid(a = sim_opts, b = chain_opts), 1, FUN = "paste", collapse = "")
names(titles) <- apply(expand.grid(a = sim_names, b = chain_names), 1, FUN = "paste", collapse = " ")

######### True vs One of the chains for each #########
chain_name = ""
# Sim 1
sim_name <- "sim1"; sim <- simulation1; df <- sim$df;
res <- readRDS(here("Simulation New Results",  paste(sim_name, "_res1", chain_name, ".RDS", sep = "")))
true <- res; true$Sj <- sim$D; true$sj_obs <- sim$d; true$rij <- sim$O
p0 <- plot_results(sim$df, true)
p1 <- plot_results(sim$df, res)

sim_name <- "sim2"; sim <- simulation2; df <- sim$df;
res <- readRDS(here("Simulation New Results",  paste(sim_name, "_res1", chain_name, ".RDS", sep = "")))
true <- res; true$Sj <- sim$D; true$sj_obs <- sim$d; true$rij <- sim$O
p2 <- plot_results(sim$df, true)
p3 <- plot_results(sim$df, res)

sim_name <- "sim3"; sim <- simulation3; df <- sim$df;
res <- readRDS(here("Simulation New Results",  paste(sim_name, "_res1", chain_name, ".RDS", sep = "")))
true <- res; true$Sj <- sim$D; true$sj_obs <- sim$d; true$rij <- sim$O
p4 <- plot_results(sim$df, true)
p5 <- plot_results(sim$df, res)



plots <- gridExtra::grid.arrange(p0, p1, p2, p3, p4, p5)
ggsave(here("Simulation Images", "Combined results plot.pdf"),
       plots,
       width = 10,
       height = 15)

######### ROC Curves - Population #########
chain_ext = "_"
par(mfrow = c(1,3))
sj_c1 <- rbind(read.csv(here("Simulation New Results", paste("Sj_sim", "sim1", "_res", 1, chain_ext, ".txt", sep = "")),header = F));
sj_c2 <- rbind(read.csv(here("Simulation New Results", paste("Sj_sim", "sim2", "_res", 1, chain_ext, ".txt", sep = "")),header = F));
sj_c3 <- rbind(read.csv(here("Simulation New Results", paste("Sj_sim", "sim3", "_res", 1, chain_ext, ".txt", sep = "")),header = F));

for(i in 2:5) {
  sj_c1 <- rbind(read.csv(here("Simulation New Results", paste("Sj_sim", "sim1", "_res", i, chain_ext, ".txt", sep = "")), header = F),sj_c1)
  sj_c2 <- rbind(read.csv(here("Simulation New Results", paste("Sj_sim",  "sim2", "_res", i, chain_ext, ".txt", sep = "")), header = F),sj_c2)
  sj_c3 <- rbind(read.csv(here("Simulation New Results", paste("Sj_sim",  "sim3", "_res", i, chain_ext, ".txt", sep = "")), header = F),sj_c3)
};

# Calculate Posterior co-clustering probability
J <- dim(sj_c1)[2]
dist_cc_mat1 <- coclusteringMatrix(sj_c1) + diag(J)
dist_cc_mat2 <- coclusteringMatrix(sj_c2) + diag(J)
dist_cc_mat3 <- coclusteringMatrix(sj_c3) + diag(J)

dist_sens_list1 <- c();dist_sens_list2 <- c();dist_sens_list3 <- c();
dist_spec_list1 <- c();dist_spec_list2 <- c();dist_spec_list3 <- c();
mcr1 <- c(); mcr2 <- c(); mcr3 <- c();
for(c in seq(0, 1.1, by = 0.01)){
  # Simulation 1
  cm <-  confusion_matrix2(simulation1, dist_cc_mat1, c, verbose = F)
  dist_sens_list1 <- c(dist_sens_list1,cm[1,1]/sum(cm[,1]))
  dist_spec_list1 <- c(dist_spec_list1, cm[2,2]/sum(cm[,2]))
  mcr1 <- c(mcr1, 1/(J*(J -1)/2)*sum(cm[1,2],cm[2,1]))

  # Simulation 2
  cm <-  confusion_matrix2(simulation2, dist_cc_mat2, c, verbose = F)
  dist_sens_list2 <- c(dist_sens_list2,cm[1,1]/sum(cm[,1]))
  dist_spec_list2 <- c(dist_spec_list2, cm[2,2]/sum(cm[,2]))
  mcr2 <- c(mcr2, 1/(J*(J -1)/2)*sum(cm[1,2],cm[2,1]))

  # Simulation 3
  cm <-  confusion_matrix2(simulation3, dist_cc_mat3, c, verbose = F)
  dist_sens_list3 <- c(dist_sens_list3,cm[1,1]/sum(cm[,1]))
  dist_spec_list3 <- c(dist_spec_list3, cm[2,2]/sum(cm[,2]))
  mcr3 <- c(mcr3, 1/(J*(J -1)/2)*sum(cm[1,2],cm[2,1]))
}
par(mfrow = c(1,1))
cols <- RColorBrewer::brewer.pal(3, "Dark2")
plot(c(0,1), c(0,1), type = "l",
     main = "Population Level ROC Curves",
     ylab = "Sensitivity",
     xlab = "1 - Specificity", col = "black", lty = 3)
lines(x = 1 - dist_spec_list1, y = dist_sens_list1, lty = 2, col = cols[1])
lines(x = 1 - dist_spec_list2, y = dist_sens_list2, lty = 2, col = cols[2])
lines(x = 1 - dist_spec_list3, y = dist_sens_list3, lty = 2, col = cols[3])

seq(0, 1.1, by = 0.01)[which(dist_spec_list1 == 1 & dist_spec_list1 == 1)]
seq(0, 1.1, by = 0.01)[which(dist_spec_list2 == 1 & dist_sens_list2 == 1)]
seq(0, 1.1, by = 0.01)[which(dist_spec_list3 == 1 & dist_spec_list3 == 1)]

######## Comparison of Epsilons ########
cols <- RColorBrewer::brewer.pal(3, "Dark2")
plot(seq(0, 1.1, by = 0.01), mcr1, type = "l",
     #main = "Effect of Epsilon on Misclassification",
     ylab = "Misclassification Rate",
     xlab = "Epsilon", lty = 1, col = cols[1])
lines(seq(0, 1.1, by = 0.01), mcr2, col = cols[2])
lines(seq(0, 1.1, by = 0.01), mcr3, col = cols[3])
legend("topright", legend = c("SD1", "SD2", 'SD3'), lty = 1, col = cols)
optimal <- seq(0, 1.1, by = 0.01)[which.min(mcr1 + mcr2 + mcr3)] #0.08
abline(v = optimal, lty = 3)


######### ROC Curves - Observational #########
sij_c1 <- rbind(read.csv(here("Simulation New Results", paste("sj_obs_sim", "sim1", "_res", 1, chain_ext, ".txt", sep = "")),header = F));
sij_c2 <- rbind(read.csv(here("Simulation New Results", paste("sj_obs_sim", "sim2", "_res", 1, chain_ext, ".txt", sep = "")),header = F));
sij_c3 <- rbind(read.csv(here("Simulation New Results", paste("sj_obs_sim", "sim3", "_res", 1, chain_ext, ".txt", sep = "")),header = F));
rij_c1 <- rbind(read.csv(here("Simulation New Results", paste("Rij_sim", "sim1", "_res", 1, chain_ext, ".txt", sep = "")),header = F));
rij_c2 <- rbind(read.csv(here("Simulation New Results", paste("Rij_sim", "sim2", "_res", 1, chain_ext, ".txt", sep = "")),header = F));
rij_c3 <- rbind(read.csv(here("Simulation New Results", paste("Rij_sim", "sim3", "_res", 1, chain_ext, ".txt", sep = "")),header = F));

for(i in 2:5) {
  sij_c1 <- rbind(read.csv(here("Simulation New Results", paste("sj_obs_sim", "sim1", "_res", i, chain_ext, ".txt", sep = "")), header = F),sij_c1)
  sij_c2 <- rbind(read.csv(here("Simulation New Results", paste("sj_obs_sim",  "sim2", "_res", i, chain_ext, ".txt", sep = "")), header = F),sij_c2)
  sij_c3 <- rbind(read.csv(here("Simulation New Results", paste("sj_obs_sim",  "sim3", "_res", i, chain_ext, ".txt", sep = "")), header = F),sij_c3)
  rij_c1 <- rbind(read.csv(here("Simulation New Results", paste("Rij_sim", "sim1", "_res", i, chain_ext, ".txt", sep = "")), header = F),rij_c1)
  rij_c2 <- rbind(read.csv(here("Simulation New Results", paste("Rij_sim",  "sim2", "_res", i, chain_ext, ".txt", sep = "")), header = F),rij_c2)
  rij_c3 <- rbind(read.csv(here("Simulation New Results", paste("Rij_sim",  "sim3", "_res", i, chain_ext, ".txt", sep = "")), header = F),rij_c3)
};

# Calculate Posterior co-clustering probability
n <- dim(sij_c1)[2]

combined1 <- matrix(paste(as.numeric(as.matrix(sij_c1)), as.numeric(as.matrix(rij_c1))),ncol = n)
combined2 <- matrix(paste(as.numeric(as.matrix(sij_c2)), as.numeric(as.matrix(rij_c2))),ncol = n)
combined3 <- matrix(paste(as.numeric(as.matrix(sij_c3)), as.numeric(as.matrix(rij_c3))),ncol = n)

obs_cc_mat1 <- coclusteringMatrix(combined1) + diag(n)
obs_cc_mat2 <- coclusteringMatrix(combined2) + diag(n)
obs_cc_mat3 <- coclusteringMatrix(combined3) + diag(n)

obs_sens_list1 <- c();obs_sens_list2 <- c();obs_sens_list3 <- c();
obs_spec_list1 <- c();obs_spec_list2 <- c();obs_spec_list3 <- c();
obs_mcr1 <- c();obs_mcr2 <- c();obs_mcr3 <- c();
for(c in seq(0, 1.1, by = 0.01)){
  # Simulation 1
  cm <-  confusion_matrix2(simulation1, obs_cc_mat1, c, verbose = F)
  obs_sens_list1 <- c(obs_sens_list1,cm[1,1]/sum(cm[,1]))
  obs_spec_list1 <- c(obs_spec_list1, cm[2,2]/sum(cm[,2]))
  obs_mcr1 <- c(obs_mcr1, 1/(n*(n -1)/2)*sum(cm[1,2],cm[2,1]))

  # Simulation 2
  cm <-  confusion_matrix2(simulation2, obs_cc_mat2, c, verbose = F)
  obs_sens_list2 <- c(obs_sens_list2,cm[1,1]/sum(cm[,1]))
  obs_spec_list2 <- c(obs_spec_list2, cm[2,2]/sum(cm[,2]))
  obs_mcr2 <- c(obs_mcr2, 1/(n*(n -1)/2)*sum(cm[1,2],cm[2,1]))

  # Simulation 3
  cm <-  confusion_matrix2(simulation3, obs_cc_mat3, c, verbose = F)
  obs_sens_list3 <- c(obs_sens_list3,cm[1,1]/sum(cm[,1]))
  obs_spec_list3 <- c(obs_spec_list3, cm[2,2]/sum(cm[,2]))
  obs_mcr3 <- c(obs_mcr3, 1/(n*(n -1)/2)*sum(cm[1,2],cm[2,1]))
}

par(mfrow = c(1,1))
cols <- RColorBrewer::brewer.pal(3, "Dark2")
plot(c(0,1), c(0,1), type = "l",
     #main = "Observation Level ROC Curves",
     ylab = "Sensitivity",
     xlab = "1 - Specificity", col = "black", lty = 3)
lines(x = 1 - obs_spec_list1, y = obs_sens_list1, lty = 1, col = cols[1])
lines(x = 1 - obs_spec_list2, y = obs_sens_list2, lty = 1, col = cols[2])
lines(x = 1 - obs_spec_list3, y = obs_sens_list3, lty = 1, col = cols[3])
legend("bottomright",
       legend = c("Random Classifier", "Simulation 1", "Simulation 2", "Simulation 3"),
       lty = c(3,1,1,1), col = c("black", cols))

seq(0, 1.1, by = 0.01)[which(obs_spec_list1 >0.9 & obs_sens_list1 >0.9)]
seq(0, 1.1, by = 0.01)[which(obs_spec_list2 >0.9 & obs_sens_list2 >0.9)]
seq(0, 1.1, by = 0.01)[which(obs_spec_list3 >0.9 & obs_sens_list3 >0.9)]

seq(0, 1.1, by = 0.01)[which(obs_mcr1 < 0.05)]
seq(0, 1.1, by = 0.01)[which(obs_mcr2 < 0.05)]
seq(0, 1.1, by = 0.01)[which(obs_mcr3 < 0.1)]



######## Comparison of Epsilons ########
cols <- RColorBrewer::brewer.pal(3, "Dark2")
plot(seq(0, 1.1, by = 0.01), obs_mcr1, type = "l",
     #main = "Effect of Epsilon on Misclassification",
     ylab = "Misclassification Rate",
     xlab = "Epsilon", lty = 1, col = cols[1])
lines(seq(0, 1.1, by = 0.01), obs_mcr2, col = cols[2])
lines(seq(0, 1.1, by = 0.01), obs_mcr3, col = cols[3])
legend("topright", legend = c("SD1", "SD2", 'SD3'), lty = 1, col = cols)
optimal <- seq(0, 1.1, by = 0.01)[which.min(obs_mcr1 + obs_mcr2 + obs_mcr3)]
abline(v = optimal, lty = 3)

plot(seq(0, 1.1, by = 0.01), dist_spec_list1, type = "l",
     main = "Effect of Epsilon on Sensitivity and Specificity",
     sub = "Population Level",
     ylab = "",
     xlab = "Epsilon", lty = 1, col = cols[1])
lines(seq(0, 1.1, by = 0.01), obs_spec_list2, col = cols[2])
lines(seq(0, 1.1, by = 0.01), obs_spec_list3, col = cols[3])
lines(seq(0, 1.1, by = 0.01), obs_sens_list1, col = cols[1], lty = 2)
lines(seq(0, 1.1, by = 0.01), obs_sens_list2, col = cols[2], lty = 2)
lines(seq(0, 1.1, by = 0.01), obs_sens_list3, col = cols[3], lty = 2)
legend(x = 0.7, y = 0.6, legend = c("Sensitivity", "Specificity"), lty = c(1,2))
legend(x = 0.7, y = 0.4, legend = c("Simulation 1", "Simulation 2", "Simulation 3"), lty = 1, col = cols)
abline(v = 0.05)

seq(0, 1.1, by = 0.01)[which.max((obs_spec_list1 + obs_sens_list1)/2 +
                                   (obs_sens_list2 + obs_spec_list2)/2 +
                                   (obs_spec_list3 +obs_sens_list3)/2)]
# 0.05
seq(0, 1.1, by = 0.01)[which.min(obs_mcr1 + obs_mcr2 + obs_mcr3)]
# 0.05
######### Heatmaps - Population #########
# Simulation 1
sim_name <- "sim1"; chain_name = ""
true_clust <- true_clusters(simulation1, "dist")
h <- hclust(as.dist(1 - dist_cc_mat1))
plot(h)
dist_idx <- order(simulation1$D)
dist_hm_df <- expand_grid(obs1 = 1:dim(dist_cc_mat1)[1], obs2 = 1:dim(dist_cc_mat1)[1])
dist_hm_df$cc_prob <- as.numeric(dist_cc_mat1)
dist_hm_df$true <- as.numeric(true_clust)

p <- create_post_cc_heatmaps(dist_hm_df, dist_idx,
                             title = paste(names(titles)[titles == paste(sim_name, chain_name, sep = "")],
                                           "Population level"))
ggsave(here("Simulation Images", paste(sim_name, chain_name, "_D_heatmaps.pdf", sep ="")),
       p, width = 10, h = 5)

# Simulation 2
sim_name <- "sim2"; chain_name = ""
true_clust <- true_clusters(simulation2, "dist")
h <- hclust(as.dist(1 - dist_cc_mat2))
plot(h)
dist_idx <- order(simulation2$D)
dist_hm_df <- expand_grid(obs1 = 1:dim(dist_cc_mat2)[1], obs2 = 1:dim(dist_cc_mat2)[1])
dist_hm_df$cc_prob <- as.numeric(dist_cc_mat2)
dist_hm_df$true <- as.numeric(true_clust)

p <- create_post_cc_heatmaps(dist_hm_df, dist_idx,
                             title = paste(names(titles)[titles == paste(sim_name, chain_name, sep = "")],
                                           "Population level"))
ggsave(here("Simulation Images", paste(sim_name, chain_name, "_D_heatmaps.pdf", sep ="")),
       p, width = 10, h = 5)

# Simulation 3
sim_name <- "sim3"; chain_name = ""
true_clust <- true_clusters(simulation3, "dist")
h <- hclust(as.dist(1 - dist_cc_mat3))
plot(h)
dist_idx <- order(simulation3$D)
dist_hm_df <- expand_grid(obs1 = 1:dim(dist_cc_mat3)[1], obs2 = 1:dim(dist_cc_mat3)[1])
dist_hm_df$cc_prob <- as.numeric(dist_cc_mat3)
dist_hm_df$true <- as.numeric(true_clust)

p <- create_post_cc_heatmaps(dist_hm_df, dist_idx,
                             title = paste(names(titles)[titles == paste(sim_name, chain_name, sep = "")],
                                           "Population level"))
ggsave(here("Simulation Images", paste(sim_name, chain_name, "_D_heatmaps.pdf", sep ="")),
       p, width = 10, h = 5)

######### Heatmaps - Observational #########
# Simulation 1
sim_name <- "sim1"; chain_name = ""
true_clust <- true_clusters(simulation1, "obs")
h <- hclust(as.dist(1 - obs_cc_mat1))
plot(h)
obs_idx <- order(simulation1$d, simulation1$O)
obs_hm_df <- expand_grid(obs1 = 1:dim(obs_cc_mat1)[1], obs2 = 1:dim(obs_cc_mat1)[1])
obs_hm_df$cc_prob <- as.numeric(obs_cc_mat1)
obs_hm_df$true <- as.numeric(true_clust)

p <- create_post_cc_heatmaps(obs_hm_df, obs_idx,
                             title = paste(names(titles)[titles == paste(sim_name, chain_name, sep = "")],
                                           "Observation level"))
ggsave(here("Simulation Images", paste(sim_name, chain_name, "_O_heatmaps.pdf", sep ="")),
       p, width = 10, h = 5)

# Simulation 2
sim_name <- "sim2"; chain_name = ""
true_clust <- true_clusters(simulation2, "obs")
h <- hclust(as.dist(1 - obs_cc_mat2))
plot(h)
obs_idx <- order(simulation2$d, simulation2$O)
obs_hm_df <- expand_grid(obs1 = 1:dim(obs_cc_mat2)[1], obs2 = 1:dim(obs_cc_mat2)[1])
obs_hm_df$cc_prob <- as.numeric(obs_cc_mat2)
obs_hm_df$true <- as.numeric(true_clust)

p <- create_post_cc_heatmaps(obs_hm_df, obs_idx,
                             title = paste(names(titles)[titles == paste(sim_name, chain_name, sep = "")],
                                           "Observation level"))
ggsave(here("Simulation Images", paste(sim_name, chain_name, "_O_heatmaps.pdf", sep ="")),
       p, width = 10, h = 5)

# Simulation 3
sim_name <- "sim3"; chain_name = ""
true_clust <- true_clusters(simulation3, "obs")
h <- hclust(as.dist(1 - obs_cc_mat3))
plot(h)
obs_idx <- order(simulation3$d, simulation3$O)
obs_hm_df <- expand_grid(obs1 = 1:dim(obs_cc_mat3)[1], obs2 = 1:dim(obs_cc_mat3)[1])
obs_hm_df$cc_prob <- as.numeric(obs_cc_mat3)
obs_hm_df$true <- as.numeric(true_clust)

p <- create_post_cc_heatmaps(obs_hm_df, obs_idx,
                             title = paste(names(titles)[titles == paste(sim_name, chain_name, sep = "")],
                                           "Observation level"))
ggsave(here("Simulation Images", paste(sim_name, chain_name, "_Oheatmaps.pdf", sep ="")),
       p, width = 10, h = 5)

######### True vs Point Estimate - need to sort label switching #########
#Simulation 1
epsilon_dist = 0.08; sim_name <- "sim1"
epsilon_obs = 0.05;
plot_actual(simulation1)
# Get results
res <- list();
for(i in 1:5)
  res[[i]] <- readRDS(here("Simulation New Results", paste(sim_name, "_res", i, chain_name, ".RDS", sep = "")))

true <- res[[1]]; true$Sj <- simulation1$D; true$sj_obs <- simulation1$d; true$rij <- simulation1$O
point_est <- res[[1]]
point_est$Sj <- cutree(hclust(as.dist(1 - dist_cc_mat1)), h = 1 - epsilon)
point_est$sj_obs <- rep(point_est$Sj, each = 10)
point_est$rij <- cutree(hclust(as.dist(1 - obs_cc_mat1)), h = 1 - epsilon)
for(i in 1:max(point_est$Sj)) {
  point_est$rij[point_est$sj_obs == i] <- exclude.empty(point_est$rij[point_est$sj_obs == i])
}
plot_results(simulation1$df, true) + plot_results(simulation1$df, point_est)

#Simulation 2
epsilon_dist = 0.08; sim_name <- "sim2"
epsilon_obs = 0.05;
# Get results
res <- list();
for(i in 1:5)
  res[[i]] <- readRDS(here("Simulation New Results", paste(sim_name, "_res", i, chain_name, ".RDS", sep = "")))
true <- res[[1]]; true$Sj <- simulation1$D; true$sj_obs <- simulation1$d; true$rij <- simulation1$O

plot_true_point_est <- function(sim, dist_cc, obs_cc, e_dist, e_obs, res){
  # Setup true
  true <- res; true$Sj <- sim$D; true$sj_obs <- sim$d; true$rij <- sim$O
  true_ordered <- label_switch(label_switch(true, "pop", sim), "obs", sim)

  # Calculate and setup point estimate
  point_est <- res
  point_est$Sj <- cutree(hclust(as.dist(1 - dist_cc)), h = 1 - e_dist)
  point_est$sj_obs <- rep(point_est$Sj, each = 10)
  point_est$rij <- cutree(hclust(as.dist(1 - obs_cc)), h = 1 - e_obs)
  for(i in 1:max(point_est$Sj)) {
    point_est$rij[point_est$sj_obs == i] <- exclude.empty(point_est$rij[point_est$sj_obs == i])
  }
  point_est_ordered <- label_switch(label_switch(point_est, "pop", sim), "obs", sim)

  # Plot
  plot_results(sim$df, true_ordered) + plot_results(sim$df, point_est_ordered)
}



plot_true_point_est(simulation3, dist_cc_mat3, obs_cc_mat3, epsilon_dist, epsilon_obs,
                    readRDS(here("Simulation New Results", paste("sim3", "_res", 1, chain_name, ".RDS", sep = ""))))
plot_true_point_est(simulation2, dist_cc_mat2, obs_cc_mat2, epsilon_dist, epsilon_obs,
                    readRDS(here("Simulation New Results", paste("sim2", "_res", 1, chain_name, ".RDS", sep = ""))))
plot_true_point_est(simulation1, dist_cc_mat1, obs_cc_mat1, epsilon_dist, epsilon_obs,
                    readRDS(here("Simulation New Results", paste("sim1", "_res", 1, chain_name, ".RDS", sep = ""))))

#Simulation 3
epsilon_dist = 0.08; sim_name <- "sim3"
epsilon_obs = 0.1;
# Get results
res <- list();
for(i in 1:5)
  res[[i]] <- readRDS(here("Simulation New Results", paste(sim_name, "_res", i, chain_name, ".RDS", sep = "")))
true <- res[[1]]; true$Sj <- simulation3$D; true$sj_obs <- simulation3$d; true$rij <- simulation3$O
point_est <- res[[1]]
point_est$Sj <- cutree(hclust(as.dist(1 - dist_cc_mat3)), h = 1 - epsilon_dist)
point_est$sj_obs <- rep(point_est$Sj, each = 10)
point_est$rij <- cutree(hclust(as.dist(1 - obs_cc_mat3)), h = 1 - epsilon_obs)
for(i in 1:max(point_est$Sj)) {
  point_est$rij[point_est$sj_obs == i] <- exclude.empty(point_est$rij[point_est$sj_obs == i])
}
plot_results(simulation3$df, true) + plot_results(simulation3$df, point_est)


######### True vs Chain Results #########
plot_results(simulation3$df, true) + plot_results(simulation3$df, point_est)


