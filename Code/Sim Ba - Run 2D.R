####### SETUP  ########
here::i_am("Code/Sim Ba - Run 2D.R")
library(here)
source(here("Code", "Multivariate mNDP Functions.R"))
source(here("Code", "MCMC result functions.R"))
source(here("Code", "Sim A - Functions.R"))

# Same for all simulations
burn_in = 2000; mcmc_iter = 1600; jumps = 5;

####### FROM SIMULATION 1 ##########
# Create small dataset (same distributions, size to fit smaller dataset closer to cell datast)
simulation1 <- simulate_data_1(30, 10, 10)
plot_actual(simulation1, 1)
#saveRDS(simulation1, here("Simulation Results", "sim1_small_dataset.RDS"))

simulation1 <- readRDS(here("Simulation Results", "sim1_small_dataset.RDS"))
df <- simulation1$df; p <- dim(df)[2] - 1

# Set parameters as per their paper, adjust nu/psi to wishart as opposed to gamma #
alpha = 1; beta = 1;
nu = 6; Psi = diag(10,p); lambda = 0.01; mu = rep(0, p)

sim1_res1 <- run_mcmc(df, p,
                      burn_in = burn_in, jumps = jumps, mcmc_iter = mcmc_iter,
                      alpha = alpha, beta = beta, S_init = 1, r_init = 1,
                      nu = nu, Psi = Psi, lambda = lambda, mu = mu, seed = NULL,
                      file_pre = "/Simulation Results/", file_post = "_sim1_res1")
saveRDS(sim1_res1, here("Simulation Results", "sim1_res1.RDS"))
plot_results(df, sim1_res1)
rm(sim1_res1)

sim1_res2 <- run_mcmc(df, p,
                      burn_in = burn_in, jumps = jumps, mcmc_iter = mcmc_iter,
                      alpha = alpha, beta = beta, S_init = 12, r_init = 2,
                      nu = nu, Psi = Psi, lambda = lambda, mu = mu, seed = NULL,
                      file_pre = "/Simulation Results/", file_post = "_sim1_res2")
saveRDS(sim1_res2, here("Simulation Results", "sim1_res2.RDS"))
plot_results(df, sim1_res2)
rm(sim1_res2)

sim1_res3 <- run_mcmc(df, p,
                      burn_in = burn_in, jumps = jumps, mcmc_iter = mcmc_iter,
                      alpha = alpha, beta = beta, S_init = 5, r_init = 4,
                      nu = nu, Psi = Psi, lambda = lambda, mu = mu, seed = NULL,
                      file_pre = "/Simulation Results/", file_post = "_sim1_res3")
saveRDS(sim1_res3, here("Simulation Results", "sim1_res3.RDS"))
plot_results(df, sim1_res3)
rm(sim1_res3)

####### FROM SIMULATION 2.R ##########
# Create small dataset (same distributions, size to fit smaller dataset closer to cell datast)
simulation2 <- simulate_data_2(30, 10, 10)
plot_actual(simulation2, 1)
#saveRDS(simulation2, here("Simulation Results", "sim2_small_dataset.RDS"))

# Only 2D version
simulation2 <- readRDS(here("Simulation Results", "sim2_small_dataset.RDS"))
df <- simulation2$df[,1:3]; p <- dim(df)[2] - 1

# Set parameters as per their paper, adjust nu/psi to wishart as opposed to gamma #
alpha = 1; beta = 1;
nu = 6; Psi = diag(10,p); lambda = 0.01; mu = rep(0, p)

sim2_res1 <- run_mcmc(df, p,
                      burn_in = burn_in, jumps = jumps, mcmc_iter = mcmc_iter,
                      alpha = alpha, beta = beta, S_init = 1, r_init = 1,
                      nu = nu, Psi = Psi, lambda = lambda, mu = mu, seed = NULL,
                      file_pre = "/Simulation Results/", file_post = "_sim2_res1")
saveRDS(sim2_res1, here("Simulation Results", "sim2_res1.RDS"))
plot_results(df, sim2_res1)
rm(sim2_res1)

sim2_res2 <- run_mcmc(df, p,
                      burn_in = burn_in, jumps = jumps, mcmc_iter = mcmc_iter,
                      alpha = alpha, beta = beta, S_init = 12, r_init = 2,
                      nu = nu, Psi = Psi, lambda = lambda, mu = mu, seed = NULL,
                      file_pre = "/Simulation Results/", file_post = "_sim2_res2")
saveRDS(sim2_res2, here("Simulation Results", "sim2_res2.RDS"))
plot_results(df, sim2_res2)
rm(sim2_res2)

sim2_res3 <- run_mcmc(df, p,
                      burn_in = burn_in, jumps = jumps, mcmc_iter = mcmc_iter,
                      alpha = alpha, beta = beta, S_init = 5, r_init = 4,
                      nu = nu, Psi = Psi, lambda = lambda, mu = mu, seed = NULL,
                      file_pre = "/Simulation Results/", file_post = "_sim2_res3")
saveRDS(sim2_res3, here("Simulation Results", "sim2_res3.RDS"))
plot_results(df, sim2_res3)
rm(sim2_res3)

####### FROM SIMULATION 3.R ##########
# Create small dataset (same distributions, size to fit smaller dataset closer to cell datast)
simulation3 <- simulate_data_3(30, 10, 10)
plot_actual(simulation3, 1)
#saveRDS(simulation3, here("Simulation Results", "sim3_small_dataset.RDS"))

# Only 2D version
simulation3 <- readRDS(here("Simulation Results", "sim3_small_dataset.RDS"))
df <- simulation3$df[,1:3]; p <- dim(df)[2] - 1

# Setup parameters
alpha = 1; beta = 1;
nu = 6; Psi = diag(10,p); lambda = 0.01; mu = rep(0, p)

sim3_res1 <- run_mcmc(df, p,
                      burn_in = burn_in, jumps = jumps, mcmc_iter = mcmc_iter,
                      alpha = alpha, beta = beta, S_init = 1, r_init = 1,
                      nu = nu, Psi = Psi, lambda = lambda, mu = mu, seed = NULL,
                      file_pre = "/Simulation Results/", file_post = "_sim3_res1")
saveRDS(sim3_res1, here("Simulation Results", "sim3_res1.RDS"))
plot_results(df, sim3_res1)
rm(sim3_res1)

sim3_res2 <- run_mcmc(df, p,
                      burn_in = burn_in, jumps = jumps, mcmc_iter = mcmc_iter,
                      alpha = alpha, beta = beta, S_init = 12, r_init = 2,
                      nu = nu, Psi = Psi, lambda = lambda, mu = mu, seed = NULL,
                      file_pre = "/Simulation Results/", file_post = "_sim3_res2")
saveRDS(sim3_res2, here("Simulation Results", "sim3_res2.RDS"))
plot_results(df, sim3_res2)
rm(sim3_res2)

sim3_res3 <- run_mcmc(df, p,
                      burn_in = burn_in, jumps = jumps, mcmc_iter = mcmc_iter,
                      alpha = alpha, beta = beta, S_init = 5, r_init = 4,
                      nu = nu, Psi = Psi, lambda = lambda, mu = mu, seed = NULL,
                      file_pre = "/Simulation Results/", file_post = "_sim3_res3")
saveRDS(sim3_res3, here("Simulation Results", "sim3_res3.RDS"))
plot_results(df, sim3_res3)
rm(sim3_res3)
