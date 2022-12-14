# Functions used when plotting and analysing results from mcmc simulation runs
library(tidyr)
library(tidyverse)
library(GGally)
library(RColorBrewer)

################################################################################
# Create a 2D plot coloured based on observational cluster &
# split by population cluster
################################################################################
plot_results <- function(dataset, res){
  p <- dim(dataset)[2] - 1
  plot.df <- as.data.frame(dataset[,2:(p + 1)])
  colnames(plot.df) <- paste("X", 1:p, sep = "")

  facet_labs <- factor(seq(1:max(res$sj_obs)), levels = seq(1:max(res$sj_obs)), ordered = T)
  names(facet_labs) <- paste("Population Cluster", facet_labs)
  S_j <- res$Sj
  sj_obs <- res$sj_obs
  rij <- res$rij
  nc = ceiling(sqrt(max(sj_obs)))
  plot.df$sj_obs = factor(sj_obs, levels = sort(unique(sj_obs)), ordered = T,
                          labels = c(paste("Population Cluster", sort(unique(sj_obs)))))
  p <- ggplot(data = plot.df, aes(x = X1, y = X2,
                                  color = factor(rij))) +
    facet_wrap(~sj_obs, ncol = nc, nrow = nc) + geom_point() + labs(color = "Observational\n Cluster",
                                                x = "", y = "") +
    theme(legend.title = element_text(size = 8),
          strip.text.x = element_text(size = 8)) +
    scale_color_brewer(palette = "Dark2")

  return(p)
}

################################################################################
# Create a 2D plot coloured based on point estimate
################################################################################
plot_point_estimate <- function(dataset, res){
  p <- dim(dataset)[2] - 1
  plot.df <- as.data.frame(dataset[,2:(p + 1)])
  colnames(plot.df) <- paste("X", 1:p, sep = "")

  facet_labs <- factor(seq(1:max(res$sj_obs)))
  names(facet_labs) <- paste("Population Cluster", facet_labs)
  S_j <- res$Sj
  sj_obs <- res$sj_obs
  rij <- res$rij

  plot.df$sj_obs = factor(sj_obs, levels = seq(1:max(sj_obs)),
                          labels = c(paste("Population Cluster", seq(1:max(sj_obs)))))
  p <- ggplot(data = plot.df, aes(x = X1, y = X2, color = factor(rij))) +
    facet_wrap(sj_obs ~.) + geom_point() + labs(color = "Observational\n Cluster",
                                                x = "", y = "") +
    theme(legend.title = element_text(size = 8),
          strip.text.x = element_text(size = 8))

  return(p)
}

################################################################################
# Create a plot of the actual simulated clusters
################################################################################
plot_actual <- function(sim, size = 0.5, d = 0){
  if(d == 0) d <- dim(sim$df)[2] - 1
  actual.df <- as.data.frame(sim$df[,2:(d + 1)], colnames(paste("V", 1:d)))
  actual.df$D <- sim$d
  actual.df$O <- sim$O
  actual.df <- unite(actual.df, "color", D:O, sep = "_", remove = F)
  if(d > 2) {
    p <- ggpairs(actual.df, columns = 1:d, aes(color = factor(color)),
                 lower = list(continuous = wrap("points", size = size)),
                 diag = list(continuous = wrap("densityDiag", alpha = 0.3)))
  }
  else { p <- ggplot(data = actual.df, aes(x = V1, y = V2, color = factor(color))) +
    facet_wrap(D ~.) + geom_point(size = size, alpha = 0.5) + labs(color = "Obs\ncluster")}
  return(p)
}

################################################################################
# Create a plot of the log probability of the data over iterations
# Color points by the number of clusters
################################################################################

plot_lp <- function(res, burnin, jump){
  lp <- res$total_l_p[(burnin + 1):length(res$total_l_p)]
  lp <- lp[seq(1:(length(lp)/jump))*jump]
  plot.df <- data.frame(it = seq(1:length(lp)),
                        lp = lp,
                        col = res$total_K[2:length(res$total_K)])
  ggplot(plot.df, aes(x = it, y = lp)) + geom_line() +
    geom_point(aes(color = factor(col)), size = 0.5) +
    labs(x = "Iteration", y = "Probability")
}

plot_lp_multiple <- function(res_list, burnin, jump){
 # base plot
  res <- res_list[[1]]
  lp <- res$total_l_p[(burnin + 1):length(res$total_l_p)]
  lp <- lp[seq(1:(length(lp)/jump))*jump]
  plot.df <- data.frame(it = seq(1:length(lp)),
                        lp = lp)
  p <- ggplot(plot.df, aes(x = it, y = lp, color = factor(1))) + geom_line() +
    labs(x = "Iteration", y = "Probability")

  for(i in 2:length(res_list)){
    res <- res_list[[1]]
    lp <- res$total_l_p[(burnin + 1):length(res$total_l_p)]
    lp <- lp[seq(1:(length(lp)/jump))*jump]
    p <- p + geom_line(aes(y = lp, color = factor(i)))
  }
  return(p)
}

################################################################################
# Map results so that output can be put into a confusion matrix for accuracy
# checking
################################################################################
map_res <- function(change_from, change_to, res){
  mapped_Sj <- res$Sj; mapped_sj_obs <- res$sj_obs
  for(i in 1:length(change_from)) {
    mapped_Sj[res$Sj == change_from[i]] = change_to[i]
    mapped_sj_obs[res$sj_obs == change_from[i]] = change_to[i]
  }
  res_mapped <- res
  res_mapped$Sj <- mapped_Sj; res_mapped$sj_obs <- mapped_sj_obs
  return(res_mapped)
}

map_vec <- function(change_from, change_to, vec){
  mapped_vec <- vec;
  for(i in 1:length(change_from)) {
    mapped_vec[vec == change_from[i]] = change_to[i]
  }
  return(mapped_vec)
}

label_switch2 <- function(res, level, sim){
  new <- res
  if(level == "pop"){
    mu <- rep(0, max(new$Sj))
    for(m in 1:max(new$Sj)) mu[m] <- mean(sim$df[new$sj_obs == m, 2])
    cat(mu, "\n")
    cat(rank(mu), "\n")
    cat(order(mu), "\n")

    new$Sj <- map_vec(1:max(new$Sj), rank(mu), new$Sj)
    new$sj_obs <- rep(new$Sj, each = 10)
    for(m in 1:max(new$Sj)) mu[m] <- mean(sim$df[new$sj_obs == m, 2])
    cat(mu, "\n")
  }

  if(level == "obs"){
    for(s in 1:max(new$Sj)){
      mu <- rep(0, max(new$rij[new$sj_obs == s]))
      for(m in 1:max(new$rij[new$sj_obs == s])) mu[m] <- mean(sim$df[new$sj_obs == s & new$rij == m, 2])
      new$rij[new$sj_obs == s] <- map_vec(1:max(new$rij[new$sj_obs == s]), rank(mu), new$rij[new$sj_obs == s] )
    }
  }
  return(new)
}

################################################################################
# Accuracy counts
################################################################################
confusion_matrix2 <- function(sim, co_clustering_matrix = NULL, epsilon = 0.01, verbose = F,
                              clusters = NULL){
  # Setup - base off provided cluster
  type <- ifelse(is.null(clusters),
                 ifelse(dim(co_clustering_matrix)[1] == length(sim$D), "dist", "obs"),
                 ifelse(length(clusters) == length(sim$D), "dist", "obs"))

  # Setup - base off co-clustering matrix
  tp = fp = fn = tn = 0
  if(is.null(clusters)){
    if(is.null(co_clustering_matrix)) {
      cat("Either cluster assignments or co_clustering_matrix must be provided.");
      return()}
    clusters <- cutree(hclust(as.dist(1 - co_clustering_matrix)), h = 1 - epsilon)}

  # For Distributional clusters
  if(type == "dist"){
    J <- length(sim$D)
    for(i in 1:(J-1)){
      for(j in (i + 1):J){
        # True Positive
        if(sim$D[i] == sim$D[j] & clusters[i] == clusters[j]){ tp <- tp + 1}
        # False Positive
        if(sim$D[i] != sim$D[j] & clusters[i] == clusters[j]){ fp <- fp + 1}
        # False Negative
        if(sim$D[i] == sim$D[j] & clusters[i] != clusters[j]){ fn <- fn + 1}
        # True Negative
        if(sim$D[i] != sim$D[j] & clusters[i] != clusters[j]){ tn <- tn + 1}
      }
    }
  }

  else if(type == "obs"){
    n <- length(sim$O)
    for(i in 1:(n-1)){
      for(j in (i + 1):n){
        # True Positive
        if(sim$d[i] == sim$d[j] & sim$O[i] == sim$O[j] & clusters[i] == clusters[j]){
          tp <- tp + 1}
        # False Positive
        if(((sim$d[i] != sim$d[j]) | (sim$O[i] != sim$O[j])) & clusters[i] == clusters[j]){
          fp <- fp + 1
        }
        # False Negative
        if(sim$d[i] == sim$d[j] & sim$O[i] == sim$O[j] & clusters[i] != clusters[j]){
          fn <- fn + 1
        }
        # True Negative
        if(((sim$d[i] != sim$d[j]) | (sim$O[i] != sim$O[j])) & clusters[i] != clusters[j]){
          tn <- tn + 1}
      }
    }
  }
  if(verbose == T) cat("TP: ", tp, "FP: ", fp, "FN: ", fn, "TN: ", tn, "\n")
  confusion_matrix <- matrix(c(tp, fn, fp, tn), 2, 2)
  return(confusion_matrix)
}


true_clusters <- function(sim, type = "obs"){
  if(type == "obs"){
    n <- dim(sim$df)[1]
    true_obs <- matrix(0, n, n)
    for(i in 1:(n)){
      for(j in i:n){
        if(sim$d[i] == sim$d[j] & sim$O[i] == sim$O[j]){
          true_obs[i,j] = true_obs[j,i] = 1
        }
      }
    }
  }

  if(type == "dist"){
    J <- length(sim$D)
    true_obs <- matrix(0, J, J)
    for(i in 1:(J)){
      for(j in i:J){
        if(sim$D[i] == sim$D[j]){
          true_obs[i,j] = true_obs[j,i] = 1
        }
      }
    }
  }

  return(true_obs)
}

######### Heatmap plots
create_post_cc_heatmaps <- function(hm_df, idx, title){
  # Plot true clusters
  p1 <- ggplot(as.data.frame(hm_df), aes(x = factor(obs1, levels =  idx, ordered = T),
                                         y = factor(obs2, levels =  idx, ordered = T), fill = true)) +
    geom_tile() + theme_bw() +
    theme(axis.ticks = element_blank(), axis.title = element_blank(),
          axis.text.x = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          legend.position = "none",
          axis.text.y = element_blank()) +
    scale_fill_gradient(low = "#ffffff", high = "#000000") +
    labs(title = title, subtitle = "True clusters")
  # Plot posterior co-clustering probability
  p2 <- ggplot(as.data.frame(hm_df), aes(x = factor(obs1, levels =  idx, ordered = T),
                                              y = factor(obs2, levels =  idx, ordered = T), fill = cc_prob)) +
    geom_tile() + theme_bw() +
    theme(axis.ticks = element_blank(), axis.title = element_blank(),
          axis.text.x = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          legend.position = "none",
          axis.text.y = element_blank()) +
    scale_fill_gradient(low = "#ffffff", high = "#000000") +
    labs(title = "", subtitle = "Posterior Co-clustering Probability", fill = "Co-clustering\n  probability",)

  return(gridExtra::grid.arrange(p1, p2, ncol = 2))
}

roc_curve2 <- function(sim, cc_mat, title = ""){
  sens_list <- c()
  spec_list <- c()
  for(c in seq(0, 1.1, by = 0.01)){
    cm <-  confusion_matrix2(sim, cc_mat, c, verbose = F)
    sens_list <- c(sens_list,cm[1,1]/sum(cm[,1]))
    spec_list <- c(spec_list, cm[2,2]/sum(cm[,2]))
  }

  cat("1 - 1: ", seq(0, 1.1, by = 0.01)[which(sens_list == 1 & spec_list == 1)])
  cat("0.8 - 0.8: ",seq(0, 1.1, by = 0.01)[which(sens_list > 0.8 & spec_list > 0.8)])

  plot(c(0,1), c(0,1), type = "l",
       main = title,
       ylab = "Sensitivity",
       xlab = "1 - Specificity", col = "red")
  lines(x = 1 - spec_list, y = sens_list, lty = 1)
  legend("bottomright", legend = c("Random Classifier",
                                   "mNDP Posterior Co-clustering"), lty = c(2, 1), col = c("red",
                                                                                           "black"))
}


