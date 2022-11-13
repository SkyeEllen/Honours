here::i_am("Code/1 Read Data.R")
library(here)
source(here("Code", "Multivariate mNDP Functions.R"))
source(here("Code", "MCMC result functions.R"))
source(here("Code", "Sim A - Functions.R"))

# Commented out to not adjust the simulated datasets
# simulation1 <- simulate_data_1(30,10,10)
# simulation2 <- simulate_data_2(30,10,10)
# simulation3 <- simulate_data_3(30,10,10)
# saveRDS(simulation1, here("Simulation New Results", "SD1.RDS"))
# saveRDS(simulation2, here("Simulation New Results", "SD2.RDS"))
# saveRDS(simulation3, here("Simulation New Results", "SD3.RDS"))
simulation1 <- readRDS(here("Simulation New Results", "SD1.RDS"))
simulation2 <- readRDS(here("Simulation New Results", "SD2.RDS"))
simulation3 <- readRDS(here("Simulation New Results", "SD3.RDS"))


chain_name = ""; burn_in = 2500; mcmc_iter = 2000; jumps = 5;
alpha = 1; beta = 1;

sim_list <- list(simulation1, simulation2, simulation3)
sim_names <- c("sim1", "sim2", "sim3")
df <- simulation2$df; p <- dim(df)[2] - 1
nu = 6; Psi = diag(10,p); lambda = 0.01; mu = rep(0, p)
S_init <- c(1,3,5,7,9); r_init <- c(1, 2, 3, 2, 1)

for(s in 1:3){
  df <- sim_list[[s]]$df; p <- dim(df)[2] - 1
  nu = 6; Psi = diag(10,p); lambda = 0.01; mu = rep(0, p)
  sim_name <- sim_names[s]
  for(i in 1:5){
    cat(paste(sim_name, "_res", i, ".RDS", sep = ""), "\n")
    res <- run_mcmc(df, p,
                    burn_in = burn_in, jumps = jumps, mcmc_iter = mcmc_iter,
                    alpha = alpha, beta = beta, S_init = S_init[i], r_init = r_init[i],
                    nu = nu, Psi = Psi, lambda = lambda, mu = mu, seed = NULL,
                    file_pre = "/Simulation New Results/", file_post = paste(sim_name, "_res", i, chain_name, sep = ""))
    saveRDS(res, here("Simulation New Results", paste(sim_name, "_res", i, chain_name, ".RDS", sep = "")))
    rm(res)
  }
}

# Increase length for simulation 1
burn_in = 2500; mcmc_iter = 3500; jumps = 5;
chain_name <- "_long"; sim_name <- "sim1"; df <- simulation1$df; p <- dim(df)[2] - 1
nu = 6; Psi = diag(10,p); lambda = 0.01; mu = rep(0, p)
for(i in 1:5){
  cat(paste(sim_name, "_res", i, chain_name, ".RDS", sep = ""), "\n")
  res <- run_mcmc(df, p,
                  burn_in = burn_in, jumps = jumps, mcmc_iter = mcmc_iter,
                  alpha = alpha, beta = beta, S_init = S_init[i], r_init = r_init[i],
                  nu = nu, Psi = Psi, lambda = lambda, mu = mu, seed = NULL,
                  file_pre = "/Simulation New Results/", file_post = paste(sim_name, "_res", i, chain_name, sep = ""))
  saveRDS(res, here("Simulation New Results", paste(sim_name, "_res", i, chain_name, ".RDS", sep = "")))
  rm(res)
}



