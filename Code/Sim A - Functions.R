library(LaplacesDemon)   # rcat
library(MASS)

####################
# Dataset based on Zuanetti et. al. 2017
# 2D dataset, but both variables draw from the same 1D distribution
# F1 =N(3,4),
# F2 = N(5,2),
# F3 =0.3N(−4,2)+0.7N(4,2),
# F4 =0.7N(−7,2) + 0.3 N(−4, 2),
# F5 = 0.4 N(−5, 2) + 0.2N(5, 2) + 0.40 N(10, 2)
# F6 = 0.15 N(−5, 2) + 0.7N(0, 2) + 0.15 N(5, 2)
####################
simulate_data_1 <- function(n_groups, I_min, I_max){
  #browser()
  ## Dataset info
  K <- 6
  J <- n_groups
  if(I_min == I_max) I_j <- rep(I_min, J)
  else I_j <- sample(seq(I_min,I_max), J, T)

  X <- matrix(0, nrow = sum(I_j), ncol = 2)

  ## Parameters for the K distributions
  sig = 0.5

  # mu
  mu_list = list(c(3),
                 c(5),
                 c(-4, 4),
                 c(-7,-4),
                 c(-5, 5, 10),
                 c(-5, 0, 5))

  # proportions
  p_list <- list(c(1),
                 c(1),
                 (c(0.3, 0.7)),
                 (c(0.7, 0.3)),
                 (c(0.4, 0.2, 0.4)),
                 (c(0.15, 0.7, 0.15)))


  D <- numeric(J)
  O <- c()
  curr <- 1
  for(j in 1:J){
    D[j] <- rcat(1, rep(1/K, K))
    for(i in 1:I_j[j]){
      O <- c(O, rcat(1, p_list[D[j]][[1]]))
      X[curr,] <- c(rnorm(1, mu_list[D[j]][[1]][O[length(O)]], sqrt(sig)),
                    rnorm(1, mu_list[D[j]][[1]][O[length(O)]], sqrt(sig)))
      curr <- curr + 1
    }

  }

  d_obs <- c(); j_obs <- c()
  for(j in 1:J) {
    d_obs <- c(d_obs, rep(D[j], I_j[j]))
    j_obs <- c(j_obs, rep(j, I_j[j]))
  }
  return(list(df = cbind(j_obs, X), O = O, D = D, d = d_obs))
}

####################
# Own Dataset with 2 well separated groups
####################
simulate_data_2 <- function(n_groups, I_min, I_max){
  K <- 2              # number of different distributions (clusters)
  J <- n_groups             # number of elements in the sample
  set.seed(1)
  if(I_min == I_max) I_j <- rep(I_min, J)
  else  I_j <- sample(c(I_min:I_max), J, T)  # number of replication for each element

  # First distribution -
  mu1a <- c(20, 20,-10)
  mu1b <- c(-5, 10, -10)
  #
  set.seed(11)
  A = matrix(runif(9), ncol = 3)
  Sigma1a <- A%*%t(A)
  set.seed(12)
  A = matrix(runif(9), ncol = 3)
  Sigma1b <-  A%*%t(A)

  mu_list1 <- list(mu1a, mu1b)
  Sigma_list1 <- list(Sigma1a, Sigma1b)
  p1 <- c(0.3, 0.7)

  # Second Distribution
  mu_list2 <- list(c(10, 0, 0), c(-10, -15, 10))
  set.seed(21)
  A = matrix(runif(9), ncol = 3)
  Sigma_2a <- A%*%t(A)
  set.seed(22)
  A = matrix(runif(9), ncol = 3)
  Sigma_2b <- round(A%*%t(A),3)
  Sigma_list2<-  list(Sigma_2a, Sigma_2b)
  p2 <- c(0.5,0.5)

  # Create List of all Mean lists and Sigma lists
  Mu_list <- list(mu_list1, mu_list2)
  Sig_list <- list(Sigma_list1, Sigma_list2)
  p_list <- list(p1, p2)


  D <- numeric(J)               # specifies which distribution each population originates from
  d <- numeric(sum(I_j))        # specifies which distribution each observation originates from
  O <- numeric(sum(I_j))        # specifies which observational distribution the observation orginates from
  X <- matrix(0, sum(I_j), 3)   # observations

  # For reproducibility
  set.seed(442)

  curr <- 1
  for(j in 1:J){
    # Get distribution
    dist <- rcat(1, rep(1/K, K))
    D[j] <- dist
    for(i in 1:I_j[j]){
      # Store
      d[curr] <- dist

      # Get observational distribution
      obsdist <- rcat(1, p_list[[dist]])
      O[curr] <- obsdist

      # Draw from distribution
      X[curr,] <- mvrnorm(1, Mu_list[[dist]][[obsdist]],
                          Sig_list[[dist]][[obsdist]])

      # Update counter
      curr <- curr + 1
    }
  }

  labs <- c()
  for(j in 1:J) labs <- c(labs, rep(j, I_j[j]))
  dataset<-matrix(c(labs,X),ncol=4) # final data set

  return(list(df = dataset, D = D, O = O, d = d))
}

####################
# Own Dataset with some groups close together
####################
simulate_data_3 <- function(n_groups, I_min, I_max){
  K <- 4              # number of different distributions (clusters)
  J <- n_groups             # number of elements in the sample
  set.seed(1)
  if(I_min == I_max) I_j <- rep(I_min, J)
  else  I_j <- sample(c(I_min:I_max), J, T)  # number of replication for each element

  # First distribution -
  mu1a <- c(19, -5,5)
  mu1b <- c(20, 5, -5)
  #
  set.seed(11)
  A = matrix(runif(9), ncol = 3)
  Sigma1a <- A%*%t(A)
  set.seed(12)
  A = matrix(runif(9, max = 3), ncol = 3)
  Sigma1b <-  A%*%t(A)

  mu_list1 <- list(mu1a, mu1b)
  Sigma_list1 <- list(Sigma1a, Sigma1b)
  p1 <- c(0.3, 0.7)

  # Second Distribution
  mu_list2 <- list(c(30, -30, 0))
  set.seed(21)
  A = matrix(runif(9), ncol = 3)
  Sigma_list2<-  list(A%*%t(A))
  p2 <- c(1)

  # Third Distribution
  mu_list3 <- list(c(-15, 8,10), c(-5, 12, 10))
  set.seed(311)
  A = matrix(runif(9, max = 3), ncol = 3)
  Sigma3a <- A%*%t(A)
  set.seed(322)
  A = matrix(runif(9), ncol = 3)
  Sigma3b <- A%*%t(A)
  Sigma_list3 <- list(Sigma3a, Sigma3b)
  p3 <- c(0.5,0.5)

  #Fourth Distribution
  mu_list4 <- list(c(-30, -10, -10), c(-27, -9, -9), c(-45, -9, -10))
  set.seed(411)
  A = matrix(rnorm(9, sd = 1), ncol = 3)
  Sigma4a <- A%*%t(A)
  set.seed(422)
  A = matrix(rnorm(9, sd = 1), ncol = 3)
  Sigma4b <- A%*%t(A)
  set.seed(433)
  A = matrix(rnorm(9, sd = 1), ncol = 3)
  Sigma4c <- A%*%t(A)
  Sigma_list4 <- list(Sigma4a, Sigma4b, Sigma4c)
  p4 <- c(0.25, 0.45, 0.3)

  # Create List of all Mean lists and Sigma lists
  Mu_list <- list(mu_list1, mu_list2, mu_list3, mu_list4)
  Sig_list <- list(Sigma_list1, Sigma_list2, Sigma_list3, Sigma_list4)
  p_list <- list(p1, p2, p3, p4)


  D <- numeric(J)               # specifies which distribution each population originates from
  d <- numeric(sum(I_j))        # specifies which distribution each observation originates from
  O <- numeric(sum(I_j))        # specifies which observational distribution the observation orginates from
  X <- matrix(0, sum(I_j), 3)   # observations

  # For reproducibility
  set.seed(442)

  curr <- 1
  for(j in 1:J){
    # Get distribution
    dist <- rcat(1, rep(1/K, K))
    D[j] <- dist
    for(i in 1:I_j[j]){
      # Store
      d[curr] <- dist

      # Get observational distribution
      obsdist <- rcat(1, p_list[[dist]])
      O[curr] <- obsdist

      # Draw from distribution
      X[curr,] <- mvrnorm(1, Mu_list[[dist]][[obsdist]],
                          Sig_list[[dist]][[obsdist]])

      # Update counter
      curr <- curr + 1
    }
  }

  labs <- c()
  for(j in 1:J) labs <- c(labs, rep(j, I_j[j]))
  dataset<-matrix(c(labs,X),ncol=4) # final data set

  return(list(df = dataset, D = D, O = O, d = d))
}

sim <- simulate_data_3(30, 10, 10)
plot_actual(sim)
