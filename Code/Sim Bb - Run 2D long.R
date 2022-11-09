####### SETUP  ########
here::i_am("Code/Sim Bb - Run 2D long.R")
library(here)
source(here("Code", "Multivariate mNDP Functions.R"))
source(here("Code", "MCMC result functions.R"))
source(here("Code", "Sim A - Functions.R"))

# Same for all simulations
burn_in = 2000; mcmc_iter = 2600; jumps = 5;

####### FROM SIMULATION 1 ##########
# Create small dataset (same distributions, size to fit smaller dataset closer to cell datast)
simulation1 <- readRDS(here("Simulation Results", "sim1_small_dataset.RDS"))
df <- simulation1$df; p <- dim(df)[2] - 1

# Set parameters as per their paper, adjust nu/psi to wishart as opposed to gamma #
alpha = 1; beta = 1;
nu = 6; Psi = diag(10,p); lambda = 0.01; mu = rep(0, p)

sim1_res1_long <- run_mcmc(df, p,
                      burn_in = burn_in, jumps = jumps, mcmc_iter = mcmc_iter,
                      alpha = alpha, beta = beta, S_init = 1, r_init = 1,
                      nu = nu, Psi = Psi, lambda = lambda, mu = mu, seed = NULL,
                      file_pre = "/Simulation Results/", file_post = "_sim1_res1_long")
saveRDS(sim1_res1_long, here("Simulation Results", "sim1_res1_long.RDS"))
plot_results(df, sim1_res1_long)
rm(sim1_res1_long)

sim1_res2_long <- run_mcmc(df, p,
                      burn_in = burn_in, jumps = jumps, mcmc_iter = mcmc_iter,
                      alpha = alpha, beta = beta, S_init = 12, r_init = 2,
                      nu = nu, Psi = Psi, lambda = lambda, mu = mu, seed = NULL,
                      file_pre = "/Simulation Results/", file_post = "_sim1_res2_long")
saveRDS(sim1_res2_long, here("Simulation Results", "sim1_res2_long.RDS"))
plot_results(df, sim1_res2_long)
rm(sim1_res2_long)

sim1_res3_long <- run_mcmc(df, p,
                      burn_in = burn_in, jumps = jumps, mcmc_iter = mcmc_iter,
                      alpha = alpha, beta = beta, S_init = 5, r_init = 4,
                      nu = nu, Psi = Psi, lambda = lambda, mu = mu, seed = NULL,
                      file_pre = "/Simulation Results/", file_post = "_sim1_res3_long")
saveRDS(sim1_res3_long, here("Simulation Results", "sim1_res3_long.RDS"))
plot_results(df, sim1_res3_long)
rm(sim1_res3_long)

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

sim2_res1_long <- run_mcmc(df, p,
                      burn_in = burn_in, jumps = jumps, mcmc_iter = mcmc_iter,
                      alpha = alpha, beta = beta, S_init = 1, r_init = 1,
                      nu = nu, Psi = Psi, lambda = lambda, mu = mu, seed = NULL,
                      file_pre = "/Simulation Results/", file_post = "_sim2_res1_long")
saveRDS(sim2_res1_long, here("Simulation Results", "sim2_res1_long.RDS"))
plot_results(df, sim2_res1_long)
rm(sim2_res1_long)

sim2_res2_long <- run_mcmc(df, p,
                      burn_in = burn_in, jumps = jumps, mcmc_iter = mcmc_iter,
                      alpha = alpha, beta = beta, S_init = 12, r_init = 2,
                      nu = nu, Psi = Psi, lambda = lambda, mu = mu, seed = NULL,
                      file_pre = "/Simulation Results/", file_post = "_sim2_res2_long")
saveRDS(sim2_res2_long, here("Simulation Results", "sim2_res2_long.RDS"))
plot_results(df, sim2_res2_long)
rm(sim2_res2_long)

sim2_res3_long <- run_mcmc(df, p,
                      burn_in = burn_in, jumps = jumps, mcmc_iter = mcmc_iter,
                      alpha = alpha, beta = beta, S_init = 5, r_init = 4,
                      nu = nu, Psi = Psi, lambda = lambda, mu = mu, seed = NULL,
                      file_pre = "/Simulation Results/", file_post = "_sim2_res3_long")
saveRDS(sim2_res3_long, here("Simulation Results", "sim2_res3_long.RDS"))
plot_results(df, sim2_res3_long)
rm(sim2_res3_long)

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

sim3_res1_long <- run_mcmc(df, p,
                      burn_in = burn_in, jumps = jumps, mcmc_iter = mcmc_iter,
                      alpha = alpha, beta = beta, S_init = 1, r_init = 1,
                      nu = nu, Psi = Psi, lambda = lambda, mu = mu, seed = NULL,
                      file_pre = "/Simulation Results/", file_post = "_sim3_res1_long")
saveRDS(sim3_res1_long, here("Simulation Results", "sim3_res1_long.RDS"))
plot_results(df, sim3_res1_long)
rm(sim3_res1_long)

sim3_res2_long <- run_mcmc(df, p,
                      burn_in = burn_in, jumps = jumps, mcmc_iter = mcmc_iter,
                      alpha = alpha, beta = beta, S_init = 12, r_init = 2,
                      nu = nu, Psi = Psi, lambda = lambda, mu = mu, seed = NULL,
                      file_pre = "/Simulation Results/", file_post = "_sim3_res2_long")
saveRDS(sim3_res2_long, here("Simulation Results", "sim3_res2_long.RDS"))
plot_results(df, sim3_res2_long)
rm(sim3_res2_long)

sim3_res3_long <- run_mcmc(df, p,
                      burn_in = burn_in, jumps = jumps, mcmc_iter = mcmc_iter,
                      alpha = alpha, beta = beta, S_init = 5, r_init = 4,
                      nu = nu, Psi = Psi, lambda = lambda, mu = mu, seed = NULL,
                      file_pre = "/Simulation Results/", file_post = "_sim3_res3_long")
saveRDS(sim3_res3_long, here("Simulation Results", "sim3_res3_long.RDS"))
plot_results(df, sim3_res3_long)
rm(sim3_res3_long)
