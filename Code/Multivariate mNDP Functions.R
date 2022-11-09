# Functions for the Multivariate mNDP
library(LaplacesDemon)   # To generate random inverse wishart
library(MASS)            # To generate random multivariate normal
library(progress)        # Progress bar for run_mcmc


################################################################################
# Sample Posterior Theta
# (Algorithm line 1/8/10)
# Equations <>
################################################################################
# Get posterior parameters for a cluster with multiple observations
get_posterior_params_cluster <- function(X, hypers){
  # Dimensions
  N <- dim(X)[1]

  # Cluster means
  Xbar <- colMeans(X)

  # Cluster covariance
  S <- t(X - Xbar)%*%(X - Xbar)

  # Posterior parameters
  mu_post <- (hypers$mu*hypers$lambda + N*Xbar)/(hypers$lambda + N)
  lambda_post <- hypers$lambda + N
  nu_post <- hypers$nu + N
  Psi_post <- hypers$Psi + S + hypers$lambda*N/lambda_post*(Xbar - hypers$mu)%*%t(Xbar - hypers$mu)

  return(list(mu = mu_post,
              lambda = lambda_post,
              nu = nu_post,
              Psi = Psi_post))
}

# Get posterior parameters for a cluster with a single observation
get_posterior_params_cluster1 <- function(X, hypers){
  X <- c(X);N <- 1
  # Posterior parameters
  mu_post <- (hypers$mu*hypers$lambda + N*X)/(hypers$lambda + N)
  lambda_post <- hypers$lambda + N
  nu_post <- hypers$nu + N
  Psi_post <- hypers$Psi + hypers$lambda*N/lambda_post*(X - hypers$mu)%*%t(X - hypers$mu)

  return(list(mu = mu_post,
              lambda = lambda_post,
              nu = nu_post,
              Psi = Psi_post))
}

# Sample from the posterior distribution for a single cluster
sample_theta_cluster <- function(X, hypers, N){
  # Get posterior parameters depending on size of cluster
  if(N > 1){
    post <- get_posterior_params_cluster(X, hypers)
  } else{
    post <- get_posterior_params_cluster1(X, hypers)
  }

  # Draw from the Inverse-Wishart distribution via wishart (Bayesian Data Analysis 2013(?))
  sig_samp <- rinvwishart(post$nu, post$Psi)
  #cat(sig_samp)
  # Ensure symmetrical (solve doesn't guarantee symmetrical even with sym input)
  if(is.symmetric.matrix(sig_samp) == F) sig_samp <- round(sig_samp, 6)

  # Draw mu from MVN
  mu_samp <- rmvn(1, post$mu, 1/post$lambda*sig_samp)
  return(list(mu= mu_samp, sig = sig_samp))
}


# Sample from the posterior distribution for all clusters
sample_posterior_theta <- function(dataset, sj, r_ij, hypers){
  # Get info from data
  N <- dim(dataset)[1]; p <- dim(dataset)[2] - 1; X <- dataset[,2:(p + 1)];

  # Identify individual clusters
  groups <- aggregate(rep(1,N), by = list(sj, r_ij), FUN = "sum")
  G <- dim(groups)[1]

  # Setup samples
  mu_samp <- matrix(0, nrow = G, ncol = p)
  Sig_samp <- array(0, c(p,p,G))

  # For each group
  for(i in 1:G){
    # Subset data
    X_kl <- X[(sj == groups$Group.1[i]) & (r_ij == groups$Group.2[i]),]

    # Get posterior draw
    samp <- sample_theta_cluster(X_kl, hypers, groups$x[i])

    # Save
    mu_samp[i,] <- samp$mu
    Sig_samp[,,i] <- samp$sig
  }

  return(list(mu = mu_samp, Sigma = Sig_samp))
}




################################################################################
# Update r_ij
# (Algorithm line 9)
# Equations <>
################################################################################
# Remove empty clusters and re-label remaining clusters accordingly
# Can be used for both distributional and observational clusters
exclude.empty <- function(clusters){
  # Empty cluster
  if(is.infinite(max(clusters))) return(clusters)

  true_max = length(table(clusters))
  need = seq(1:true_max)
  while(true_max < max(clusters)){
    labels = as.numeric(as.character(data.frame(table(clusters))[,1]))
    dif <- which(need != labels)
    clusters[clusters > dif[1]] = clusters[clusters > dif[1]] - 1
  }
  return(clusters)
}

# Update the observational cluster for a given group j
update_observational_cluster <- function(dataset, Sj, Rij, sj, Lk, theta.samp, beta, hypers, j){
  # Current groups (j,i)
  groups = aggregate(dataset, by = list(sj, Rij), FUN = "sum")[,1:2]

  # Distributional cluster
  k <- Sj[j]

  # Indices of observations in kth distributional cluster
  indices <- which(sj == k)

  # Observations in the kth distributional cluster
  Xk <- dataset[indices, 2:dim(dataset)[2]]

  # Observational cluster for kth distributional cluster
  Rk <- Rij[indices]

  # Counter
  curr <- 1

  # Calculate probability of rij = l for each observation
  for(i in indices){
    # Remove current observation from set
    Rk.minus = Rk[-curr]

    # Count for each observational cluster
    # Consider the case that removing the observation results in empty cluster
    mlk.minus = get_mlk(Lk[k], Rk.minus)

    prob_rij <- get_prob_rij(Xk[curr,], Lk[k], k,
                             groups, theta.samp, mlk.minus, hypers, beta)

    # Reassign observational cluster
    Rij[i]<-rcat(1, prob_rij)

    # Update observational cluster within kth distributional cluster
    Rk<-Rij[indices]

    # Increment counter
    curr<-curr+1
  }

  # Exclude empty observational clusters
  Rij[indices] <- Rk
  Rij <- exclude.empty(Rij)

  return(Rij)
}

################################################################################
# Calculate probability of r_ij
# Part of calculating acceptance probability
################################################################################
# Line 1 - probability based on theta sample
pr_ij_mvn <- function(obs, m, mu, Sigma){
  m*dmvn(obs, mu, Sigma)
}

# Line 2 - drawn from the multivariate t distribution
pr_ij_t <- function(obs, b, hypers){
  b*dmvt(obs,
         mu = hypers$mu,
         S = (hypers$lambda + 1)/(hypers$lambda*(hypers$nu - length(obs) + 1))*hypers$Psi,
         df = hypers$nu - length(obs) + 1)
}

# Calculate the probability of r_ij for a given set of observations
get_prob_rij <- function(obs, Lk, k, groups, theta.samp, mlk, hypers, beta){
  # Initalise
  prob_rij <- numeric(Lk)

  # Get probability for each cluster l
  for(l in 1:Lk){
    # Location of parameters w/in parameter samples
    idx = which(groups$Group.1 == k & groups$Group.2 == l)

    if(length(idx) != 0){
      # Parameters
      params = list(mu = theta.samp$mu[idx,], Sigma = theta.samp$Sigma[,,idx])
      # Probability
      prob_rij[l] = pr_ij_mvn(obs, mlk[l], params$mu, params$Sigma)
    }
    else {
      # Probability for empty group
      prob_rij[l] = pr_ij_t(obs, beta, hypers)
    }
  }

  # Probability of a new cluster
  prob_rij <- c(prob_rij, pr_ij_t(obs, beta, hypers))

  # Scale probabilities to one
  prob_rij = prob_rij/sum(prob_rij)

  return(prob_rij)
}




################################################################################
# Build candidate move
################################################################################
# Calculate mlk for a given k
# given the number of observational clusters Lk and observational cluster indicators
get_mlk <- function(Lk, Rk){
  mlk <- numeric(Lk)
  for(l in 1:Lk){ mlk[l] = length(Rk[Rk == l]) }
  return(mlk)
}

# Build proposal observational cluster allocation if
# proposal distributional cluster is new
build_obs_cluster_Knew <- function(data, hypers_DP, hypers_G){
  n <- dim(data)[1]; p <- dim(data)[2]
  l_prob_rj_prop <- l_q_prop <- 0
  # First observation is necessarily new observational cluster
  Rij_prop <- numeric(n); Rij_prop[1] <- 1

  # Remaining observations
  for(i in 2:n){
    if(n == 1) break
    # Number of observations in each observational cluster
    # (excluding observations that haven't been reassigned)
    mlk <- c(table(Rij_prop[Rij_prop != 0]), hypers_DP$beta)

    # Calculate probability for each cluster and scale to 1
    prob_rij <- mlk*pr_ij_t(data[i,], 1, hypers_G)
    prob_rij <- prob_rij/sum(prob_rij)

    # Draw observational cluster
    Rij_prop[i] <- rcat(1, prob_rij)

    # Scale mlk
    mlk <- mlk/sum(mlk)

    # Marginal probability of r
    l_prob_rj_prop <- l_prob_rj_prop + log(mlk[Rij_prop[i]])

    # Transition function
    l_q_prop <- l_q_prop + log(prob_rij[Rij_prop[i]])
  }

  return(list(r_ij = Rij_prop, l_prob_rj_prop=l_prob_rj_prop, l_q_prop=l_q_prop))
}

# Build proposal observational cluster allocation if
# proposal distributional cluster is old
build_obs_cluster_Kold <- function(Xj, j, Ij, sj_obs, r_ij, Sprop, theta_samp, groups,
                                   hypers_G, hypers_DP){
  n <- dim(Xj)[1];
  rk <- r_ij[sj_obs == Sprop]; Lk <- max(rk)
  mlk <- c(get_mlk(Lk, rk), hypers_DP$beta)

  rj_prop <- numeric(n); l_prob_rj_prop <- 0; l_q_prop <- 0;
  for(i in 1:Ij[j]){
    # Get probabilities
    prob_rij <- get_prob_rij(Xj[i,], Lk, Sprop, groups, theta_samp, mlk[1:Lk], hypers_G, hypers_DP$beta)

    rj_prop[i] <- rcat(1, prob_rij)
    # Marginal probability of r
    l_prob_rj_prop <- l_prob_rj_prop + log(mlk[rj_prop[i]]/sum(mlk))

    # Transition function
    l_q_prop <- l_q_prop + log(prob_rij[rj_prop[i]])
  }

  return(list(r_ij = rj_prop, l_prob_rj_prop=l_prob_rj_prop, l_q_prop=l_q_prop))
}

# Build proposal observational cluster allocation if
# proposal distributional cluster is new
build_proposal <- function(dataset, Sj, sj_obs, j, K, hypers_DP, Ij, hypers_G, theta_samp, r_ij){
  ##### Setup for calculating probabilities
  # Current distributional cluster
  S_curr <- Sj[j]
  # Current observational clusters (built alongside proposal)
  j_indices <- seq(cumsum(Ij)[j] - Ij[j] + 1, cumsum(Ij)[j])
  rj <- r_ij[j_indices]
  r_ij_curr <- numeric(Ij[j])

  # log probabilities (sum)
  l_prob_rj_prop <- l_q_prop <- 0

  # Observations in jth cluster
  Xj <- dataset[j_indices, 2:dim(dataset)[2]]

  ## Proposal distributional cluster ##
  # Number of groups in each distributional cluster (excluding jth group)
  nk_minj <- numeric(K)
  for(k in 1:K) nk_minj[k] <- sum(Sj[-j] == k)

  # Include probability of new distributional cluster
  nk_minj <- c(nk_minj, hypers_DP$alpha)

  # Force it to be a different cluster
  nk_minj[Sj[j]] <- 0

  # Scale to become probability for each cluster
  prob_k <- nk_minj/sum(nk_minj)

  # Get proposal cluster
  Sprop <- rcat(1, prob_k)

  ## Proposal observational clusters
  # If Sprop is in a current distributional cluster
  if(Sprop <= K){
    groups = aggregate(dataset, by = list(sj_obs, r_ij), FUN = "sum")[,1:2]
    obs_clust <- build_obs_cluster_Kold(Xj, j, Ij, sj_obs, r_ij, Sprop, theta_samp, groups,
                                        hypers_G, hypers_DP)
  }

  # If Sprop is a new distributional cluster (K + 1)
  if (Sprop > K) {
    obs_clust <- build_obs_cluster_Knew(Xj, hypers_DP, hypers_G)
  }

  # Combine for full candidate
  Sj_prop <- Sj
  Sj_prop[j] <- Sprop

  # At observational level
  sj_obs_prop <- sj_obs
  sj_obs_prop[j_indices] <- rep.int(Sprop, Ij[j])

  # observational clusters
  rij_prop <- r_ij
  rij_prop[j_indices] <- obs_clust$r_ij

  proposal <- list(Sj = Sj_prop,
                   sj_obs = sj_obs_prop,
                   rij = rij_prop,
                   l_prob_rj_prop = obs_clust$l_prob_rj_prop,
                   l_q_prop=obs_clust$l_q_prop )

  return(proposal)
}


################################################################################
# Calculate log probability of the current allocation for group j
################################################################################
log_prob_curr <- function(dataset, j, Sj, sj_obs, Ij, Rij, K, Lk, theta.samp, hypers_DP, hypers_G){
  ##### Setup for calculating probabilities
  # Current observational clusters (built alongside proposal)
  j_indices <- seq(cumsum(Ij)[j] - Ij[j] + 1, cumsum(Ij)[j])
  Rj <- Rij[j_indices]
  Rij_curr <- numeric(Ij[j])

  # log probabilities (sum)
  l_prob_rj_curr <- 0; l_q_curr <- 0

  # Observations in jth cluster
  Xj <- dataset[j_indices, 2:dim(dataset)[2]]

  # Number of groups in each distributional cluster (excluding jth group)
  nk_minj <- numeric(K)
  for(k in 1:K) nk_minj[k] <- sum(Sj[-j] == k)

  # Number of groups remaining in jth group's cluster without jth group
  nj <- nk_minj[Sj[j]]

  # Treat as new cluster if the current cluster is empty
  if(nj == 0){
    # First (should) be 1
    Rij_curr[1] <- Rj[1]
    l_prob_rj_curr = l_prob_rj_curr
    l_q_curr = l_q_curr

    # Remaining observations
    for(i in 2:Ij[j]){
      # Number of observations in each observational cluster
      # (excluding observations that haven't been reassigned)
      mlk <- get_mlk(Lk[Sj[j]], Rij_curr)

      # Chance of new cluster
      mlk <- c(mlk, hypers_DP$beta)

      # Calculate probability for each cluster and scale to 1
      prob_rij <- mlk*pr_ij_t(Xj[i,], 1, hypers_G)
      prob_rij <- prob_rij/sum(prob_rij)

      # Scale mlk
      mlk <- mlk/sum(mlk)

      # Get current observational cluster
      Rij_curr[i] <- Rj[i]

      # Marginal probability of r
      # New cluster
      if(mlk[Rij_curr[i]] == 0){
        l_prob_rj_curr <- l_prob_rj_curr + log(mlk[length(mlk)])
        l_q_curr <- l_q_curr + log(prob_rij[length(mlk)])
      } else {
        l_prob_rj_curr <- l_prob_rj_curr + log(mlk[Rij_curr[i]])
        l_q_curr <- l_q_curr + log(prob_rij[Rij_curr[i]])
      }
    }
  }

  # Get probability based on remaining clusters if not empty
  if(nj > 0){
    # Indicator variables without jth cluster
    j_indices <- seq(cumsum(Ij)[j] - Ij[j] + 1, cumsum(Ij)[j])
    R_minusj <- Rij[-j_indices]
    sj_obs_minusj <- sj_obs[-j_indices]

    # Subset to observations in kth distributional cluster
    k_indices <- which(sj_obs_minusj == Sj[j])
    Rk_minusj <- R_minusj[k_indices]
    Lk_minus <- max(Rk_minusj)

    # Current groups (j,i)
    groups = aggregate(dataset, by = list(sj_obs, Rij), FUN = "sum")[,1:2]

    for(i in 1:Ij[j]){
      ## Calculate probabilities
      curr <- setup(Xj[i,],Rk_minusj, Rij_curr, Lk_minus, Sj[j],
                     groups, theta.samp, hypers_G, hypers_DP$beta)

      # Get current
      Rij_curr[i] <- Rj[i]
      # New cluster
      if(Rij_curr[i] > (Lk_minus + 1)){
        l_prob_rj_curr <- l_prob_rj_curr + log(curr$mlk[length(curr$mlk)])
        l_q_curr <- l_q_curr + log(curr$prob_rij[length(curr$prob_rij)])
      }
      # Empty cluster (same as new cluster)
      else if(curr$mlk[Rij_curr[i]] == 0){
        l_prob_rj_curr <- l_prob_rj_curr + log(curr$mlk[length(curr$mlk)])
        l_q_curr <- l_q_curr + log(curr$prob_rij[length(curr$prob_rij)])
      }
      # Old cluster
      else {
        l_prob_rj_curr <- l_prob_rj_curr + log(curr$mlk[Rij_curr[i]])
        l_q_curr <- l_q_curr + log(curr$prob_rij[Rij_curr[i]])
      }
    }
  }

  return(list(l_prob_rj_curr = l_prob_rj_curr,
              l_q_curr = l_q_curr))
}

# Setup the probabilities for a given cluster k with jth group removed
setup <- function(obs, Rk_minusj, Rj, Lk_minus, k, groups, theta_samp, hypers_G, beta){
  # Add on any redrawn observations
  Rk <- c(Rk_minusj, Rj[Rj > 0])

  # Count number of observations in each cluster
  mlk <- get_mlk(max(Rk), Rk)

  # Calculate probability for each
  prob_rij <- numeric(max(Rk))
  for(l in 1:Lk_minus){
    # Location of parameters w/in parameter samples
    idx = which(groups$Group.1 == k & groups$Group.2 == l)

    if(length(idx) != 0){
      # Probability
      prob_rij[l] = mlk[l]*dmvn(obs, theta_samp$mu[idx,],
                                theta_samp$Sigma[,,idx], log = F)
    }
    else {
      # Probability for empty group
      prob_rij[l] = 0
    }
  }

  if(max(Rk) > Lk_minus){
    for(l in (Lk_minus + 1):max(Rk)){
      prob_rij[l] = pr_ij_t(obs, mlk[l], hypers_G)
    }
  }

  prob_rij <- c(prob_rij, pr_ij_t(obs, beta, hypers_G))
  prob_rij <- prob_rij/sum(prob_rij)

  mlk <- c(mlk, beta)/sum(mlk, beta)
  return(list(mlk = mlk, prob_rij = prob_rij))
}

##############################################################################
# Calculate the log marginal likelihood of X_ij for a given j
#
# Input:
# Xjl - data within the wanted subset
# hypers - posterior theta parameters calculated without the jth group
# Output: log density
#
##############################################################################

dmarginal_likelihood_j <- function(Xj, groups, S, rj, post_minusj, hypers_G){
  lmarg_like <- 0
  for(l in unique(rj)){
    l_indices <- which(rj == l)

    param_idx <- which(groups$Group.1 == S & groups$Group.2 == l)

    # Use parameters from posterior
    if(length(param_idx == 1)){
      hypers <- list(mu = post_minusj$mu[param_idx,],
                     nu = post_minusj$nu[param_idx],
                     lambda = post_minusj$lambda[param_idx],
                     Psi = post_minusj$Psi[,,param_idx])
      new <- dmarginal_likelihood(Xj[l_indices,], hypers)
      if(length(new) == 0) {
        print(param_idx)
      }
      lmarg_like = lmarg_like + dmarginal_likelihood(Xj[l_indices,], hypers)
    } else {
      # Use parameters from prior if group doesn't exist in -j
      lmarg_like = lmarg_like + dmarginal_likelihood(Xj[l_indices,], hypers_G)
    }
  }

  return(lmarg_like)
}

##############################################################################
# Calculate the log marginal likelihood of X_ij for a given j and i = l
#
# Input:
# Xjl - data within the wanted subset
# hypers - posterior theta parameters calculated without the jth group
# Output: log density
#
##############################################################################
dmarginal_likelihood<-function(Xjl, hypers){
  p <- dim(Xjl)[2]
  m_lkj <- dim(Xjl)[1]

  if(is.null(m_lkj)){
    m_lkj <- 1
    p <- length(Xjl)
    mu_lkj <- Xjl
    Sigma_lkj <- matrix(0,p,p)
  }
  else{
    mu_lkj <- colMeans(Xjl)
    Sigma_lkj <-  t(Xjl - mu_lkj)%*%(Xjl - mu_lkj)
  }

  nu_j <- hypers$nu + m_lkj
  lambda_j <- hypers$lambda + m_lkj
  mu_diff <- mu_lkj - hypers$mu
  mu_diff <- as.matrix(as.numeric(mu_diff))
  Psi_j <- hypers$Psi + Sigma_lkj + (hypers$lambda*m_lkj)/lambda_j*((mu_diff)%*%t(mu_diff))

  marg_likelihood <- -m_lkj/2*log(2*pi) + m_lkj*p/2*log(2) +
    lgamma(nu_j/2) - lgamma(hypers$nu/2) +
    p/2*log((hypers$lambda)/(lambda_j)) +
    hypers$nu/2*log(det(hypers$Psi)) - nu_j/2*log(det(Psi_j))

  return(marg_likelihood)
}

#####
# Calculate the probability of accepting candidate
# Input parameters:
  # dataset: dataset including group labels
  # Sj: distributional cluster labels at group levels
  # sj_obs: distributional cluster labels at observational levels
  # r_ij: observational cluster labels
  # candidate: list of candidate values
  # curr_probs: probabilities related to current cluster allocation
  # Ij: number of observations in each group
  # j: group being updated
  # hypers_G: hyperparameters of the base distribution
#####
# Get posterior parameters for each cluster (internal function)
get_posterior_params <- function(X, sj_obs, r_ij, hypers, groups){
  G <- dim(groups)[1]; p <- dim(X)[2];
  mu_list <- matrix(0, nrow = G, ncol = p)
  Psi_list <- array(0, c(p,p,G))
  nu_list <- numeric(G)
  lambda_list <- numeric(G)

  # For each group
  for(i in 1:G){
    # Subset data
    X_kl <- matrix(X[(sj_obs == groups$Group.1[i]) & (r_ij == groups$Group.2[i]),], ncol = p)

    if(dim(X_kl)[1] > 1){
      post <- get_posterior_params_cluster(X_kl, hypers)
    } else{
      post <- get_posterior_params_cluster1(X_kl, hypers)
    }

    # Save
    mu_list[i,] <- post$mu; Psi_list[,,i] <- post$Psi; nu_list[i] <- post$nu; lambda_list[i] <- post$lambda
  }

  return(list(mu = mu_list, Psi = Psi_list, nu = nu_list, lambda = lambda_list))
}

# Calculate acceptance probability
prob_accept <- function(dataset, Sj, sj_obs, r_ij, candidate, curr_probs, Ij, j,
                        hypers_G){
  # Subset data
  j_indices <- seq(cumsum(Ij)[j] - Ij[j] + 1, cumsum(Ij)[j])
  Xj <- dataset[j_indices, 2:dim(dataset)[2]]

  # Cluster indicators current
  S <- Sj[j]; s <- sj_obs[j_indices]; rj <- r_ij[j_indices]

  # Cluster indicators candidate
  S_prop <- candidate$Sj[j]; s_prop <- candidate$sj_obs[j_indices];
  rj_prop <- candidate$rij[j_indices]

  # Posterior without jth group
  X_minusj <- dataset[-j_indices, 2:dim(dataset)[2]]
  groups <- aggregate(rep(1, nrow(X_minusj)), by = list(sj_obs[-j_indices], r_ij[-j_indices]), FUN = "sum")[,1:2]
  post_minusj <- get_posterior_params(X_minusj, sj_obs[-j_indices],
                                      r_ij[-j_indices], hypers_G, groups)

  # Log marginal likelihood of Xj
  lmarg_curr <-  dmarginal_likelihood_j(Xj, groups, S, rj, post_minusj, hypers_G)
  lmarg_prop <-  dmarginal_likelihood_j(Xj, groups, S_prop, rj_prop, post_minusj, hypers_G)

  # Combined with other probabilities
  p_accept <- exp(candidate$l_prob_rj_prop + lmarg_prop + curr_probs$l_q_curr -
                    (curr_probs$l_prob_rj_curr + lmarg_curr + candidate$l_q_prop))
  return(list(p_accept = p_accept, lmarg_prop = lmarg_prop, lmarg_curr = lmarg_curr))
}

##############################################################################
# Calculate the log marginal likelihood of all data (for convergence testing)
# Input:
# dataset -
# hypers - posterior theta parameters calculated without the jth group
# Output: log density
#
##############################################################################
y_log_prob <- function(dataset, r_ij, sj_obs, theta_samp, hypers_G){
  #browser()
  N <- dim(dataset)[1]; p <- dim(dataset)[2]
  groups <- aggregate(rep(1,N), by = list(sj_obs, r_ij), FUN = "sum")[,1:2]
  l_prob <- 0
  for(i in 1:N){
    param_idx <- which(groups$Group.1 == sj_obs[i] & groups$Group.2 == r_ij[i])
    l_prob = l_prob + dmvn(dataset[i, 2:p], theta_samp$mu[param_idx,],
                           theta_samp$Sigma[,,param_idx], log = T)
  }
  return(l_prob)
}
#########################
# Run MCMC
# dataset = data, p = number of variables in data,
# S_init = number of distributional clusters to start with (integer)
# r_init = number of observational clusters to start with (integer)
# alpha, beta = hyperparameters for DP
# mu, lambda, nu, Psi = hyperparameters for base distribution (Normal Inverse Wishart)
# burn_in = mcmc burn in iterations
# jumps = length of jump
# file_extension - for output files of mcmc iteration info
#########################
run_mcmc <- function(dataset, p, S_init = 1, r_init = 1,
                     alpha = 1, beta = 1,
                     mu = rep(0, p), lambda = 0.01, nu = 10, Psi = diag(p),
                     burn_in = 1000, jumps = 10, mcmc_iter = 1000,
                     file_pre = "/", file_post = "", debug = FALSE,
                     progress = TRUE, output = TRUE, seed = 1234){

  #########################
  # Model initialization
  #########################
  wd <- getwd()

  # Descriptive
  J <- max(dataset[,1])                           # No. of groups
  I_j <- numeric(J)                               # No. of obs in each group
  for (i in 1:J) I_j[i]<-sum(dataset[,1]==i)
  n_obs <- nrow(dataset)                          # Total no. of obs
  p <- ncol(dataset) - 1                          # Total no. of predictors

  set.seed(NULL)
  S_j <- exclude.empty(sample(1:S_init, J, T))      # Distributional cluster indicator (group level)
  sj_obs <- numeric()         # Distributional cluster indicator (observational level)
  for (j in 1:length(S_j))    sj_obs <- c(sj_obs, rep(S_j[j], I_j[j]))

  r_ij <- rep(1, n_obs)       # Observational cluster indicator
  if(S_init > 1){
    for(i in 1:S_init){
      indx <- which(sj_obs == 1)
      r_ij[indx] <- exclude.empty(sample(1:r_init, length(indx), T))
    }
  }

  K <- length(table(S_j))     # Number of observed distributional clusters
  L_k <- numeric(K)           # Number of observational clusters per dc
  for (k in 1:K)              L_k[k] <- max(r_ij[which(sj_obs == k)])

  # Hyperparameters - DP
  hypers_DP <- list(alpha = alpha,           # concentration parameter of outer DP
                    beta = beta)            # concentration parameter of inner DP

  # Hyperparameters - Normal Inverse Wishart Distribution
  hypers_G <- list(mu = mu,       # mean vector for Normal
                   lambda = lambda,           # scaling parameter for cov matrix
                   nu = nu,              # degrees of freedom for Inv Wishart
                   Psi = Psi)     # covariance matrix for Inv Wishart

  #########################
  # MCMC initialization
  #########################
  # Initialise sample of thetas
  # Set seed for replication
  set.seed(seed)
  theta_samp <- sample_posterior_theta(dataset, sj_obs, r_ij, hypers_G)

  total_p_accept <- 1         # Vector of probability of acceptance at each iteration
  total_accepted <- 0         # Vector: 0 if accepted, 1 if not
  total_a <- c(); total_j <- c(); total_new_Sj <- c(); total_K <- c(K); log_p_total <- c()

  # MCMC iteration values
  total_iter <- burn_in + mcmc_iter*jumps

  # Seed for reproducibility
  set.seed(seed)

  #########################
  # Run MCMC
  #########################
  if(progress == TRUE) {
    pb <- progress_bar$new(format = "[:bar] :percent || Estimated time remaining: :eta",
                           total = total_iter, complete = "-", incomplete = " ",
                           current = ">", clear = FALSE, width = 100)
  }


  for(it in 1:total_iter){
    # Select j to update
    j <- sample(1:J, 1)
    total_j <- c(total_j, j)

    # Obtain candidate move for distributional cluster
    cand <- build_proposal(dataset, S_j, sj_obs, j, K, hypers_DP,
                           I_j, hypers_G, theta_samp, r_ij)
    if(debug == T) browser()
    #plot_cand(dataset, sj_obs, S_j, j, cand)
    curr_p <- log_prob_curr(dataset, j, S_j, sj_obs, I_j, r_ij, K, L_k, theta_samp, hypers_DP,
                            hypers_G)

    # Calculate acceptance probability & store
    accept_probs <- prob_accept(dataset, S_j, sj_obs, r_ij, cand, curr_p, I_j, j, hypers_G)
    p_accept <- accept_probs$p_accept
    total_p_accept <- c(total_p_accept, p_accept)

    # Calculate comparison value
    a <- runif(1)
    total_a <- c(total_a, a)

    # Accept move
    if (a < p_accept){
      if(output == TRUE) cat("\n", j, K)

      total_accepted <- c(total_accepted, 1)
      total_new_Sj <- c(total_new_Sj, cand$Sj[j])

      # Update indicators with candidate values
      S_j <- cand$Sj
      sj_obs <- cand$sj_obs # ( will be updated below anyway...?)
      r_ij <- cand$rij

      # Remove empty distributional clusters
      S_j <- exclude.empty(S_j)
      sj_obs <- numeric()
      for(i in 1:length(S_j)) sj_obs <- c(sj_obs, rep(S_j[i], I_j[i]))

      # Remove empty observational clusters
      K <- max(S_j)
      L_k <- numeric(K)
      for(k in 1:K){
        obs <- which(sj_obs == k)
        r_ij[obs] <- exclude.empty(r_ij[obs])
        L_k[k] <- max(r_ij[obs])
      }

      # Update parameters
      theta_samp <- sample_posterior_theta(dataset, sj_obs, r_ij, hypers_G)

      # Update observational clusters
      r_ij <- update_observational_cluster(dataset, S_j, r_ij, sj_obs, L_k, theta_samp,
                                           hypers_DP$beta, hypers_G, j)


      # Update L_k observational clusters in each distributional cluster
      for(k in 1:K){
        obs <- which(sj_obs == k)
        r_ij[obs] <- exclude.empty(r_ij[obs])
        L_k[k] <- max(r_ij[obs])
      }

      if(output == TRUE) cat("\n", table(S_j))

      # Sample parameters under new observational clusters
      theta_samp <- sample_posterior_theta(dataset, sj_obs, r_ij, hypers_G)

    }

    if(a >= p_accept) {
      total_accepted <- c(total_accepted, 0)
    }

    # Log p(y | r, s, theta)
    log_p_total <- c(log_p_total, y_log_prob(dataset, sj_obs, r_ij, theta_samp, hypers_G))
    total_K <- c(total_K, K)

    # Save parameters intermittently for chains
    if(it > burn_in & it%%jumps == 0){
      if(it >= burn_in + jumps) append = T
      # Distributional clusters
      write(sj_obs,  ncolumns = length(sj_obs), file=paste(wd, file_pre, "sj_obs_sim", file_post,".txt",sep=""), sep = ",", append=T)
      write(S_j,  ncolumns = length(S_j), file=paste(wd, file_pre, "Sj_sim", file_post,".txt",sep=""), sep = ",", append=T)
      write(K,  ncolumns = 1, file=paste(wd,file_pre,"K_sim", file_post,".txt",sep=""),sep = ",", append=T)

      # Observational clusters
      write(r_ij,  ncolumns = length(r_ij), file=paste(wd,file_pre,"Rij_sim", file_post,".txt",sep=""),sep = ",", append=T)
      write(L_k,  ncolumns = length(L_k), file=paste(wd,file_pre,"Lk_sim", file_post,".txt",sep=""),sep = ",", append=T)

      # Parameters
      write(theta_samp$mu, ncolumns = length(theta_samp$mu), file=paste(wd,file_pre,"mu_sim", file_post,".txt",sep=""),sep = ",", append=T)
      write(theta_samp$Sigma, ncolumns = length(theta_samp$Sigma), file=paste(wd,file_pre,"Sigma_sim", file_post,".txt",sep=""),sep = ",", append=T)
    }
    if(progress == TRUE) pb$tick()

  }

  return(list(Sj = S_j, rij = r_ij, total_j = total_j, sj_obs = sj_obs,
              total_p_accept = total_p_accept, total_accepted = total_accepted,
              total_K = total_K, total_l_p = log_p_total))
}


