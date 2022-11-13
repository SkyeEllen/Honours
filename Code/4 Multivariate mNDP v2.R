# Application of the mNDP to both datasets & visualisations/analysis of results
here::i_am("Code/4 Multivariate mNDP.R")
library(here)
library(tidyverse)
source(here("Code", "Multivariate mNDP Functions.R"))
source(here("Code", "MCMC result functions.R"))

cell_umap <- readRDS(here("RNA Splicing Data", "Cell UMAP.RDS"))
tissue_umap <- readRDS(here("RNA Splicing Data", "Tissue UMAP.RDS"))
bio_source_all <- readRDS(here("RNA Splicing Data", "Bio_source_all.RDS"))

cell_source <-  readRDS(here("RNA Splicing Data", "Cell source.RDS"))
tissue_source <-  readRDS(here("RNA Splicing Data", "Tissue source.RDS"))

# Set up df for mNDP
cell_source$Bio_source_num <- as.numeric(as.factor(cell_source$Biological_source))
cell_df <- left_join(cell_source,
                     data.frame(cell_umap, RNA_number_id = rownames(cell_umap)))
rownames(cell_df) <- cell_df$RNA_number_id
cell_df$RNA_number_id <- NULL; cell_df$Biological_source <- NULL;
cell_df <- cell_df[order(cell_df$Bio_source_num),]
saveRDS(cell_df, here("RNA Splicing Data", "cell group df.RDS"))

# Set up tissue df for mNDP
tissue_source$Cell_sys_num <- as.numeric(as.factor(tissue_source$Cell_system))
tissue_df <- left_join(tissue_source,
                       data.frame(tissue_umap, RNA_number_id = rownames(tissue_umap)))
rownames(tissue_df) <- tissue_df$RNA_number_id
tissue_df$RNA_number_id <- NULL; tissue_df$Biological_source <- NULL; tissue_df$Cell_system <- NULL
tissue_df <- tissue_df[!(tissue_df$Cell_sys_num %in% which(table(tissue_df$Cell_sys_num) == 1)),]
tissue_df$Cell_sys_num <- exclude.empty(tissue_df$Cell_sys_num)
tissue_df <- tissue_df[order(tissue_df$Cell_sys_num),]
saveRDS(tissue_df, here("RNA Splicing Data", "tissue group df.RDS"))

####### cell df ##########
df <- as.matrix(cell_df); p <- dim(df)[2] - 1

# Set parameters as per their paper, adjust nu/psi to wishart as opposed to gamma #
alpha = 0.5; beta = 2; burn_in = 2000; mcmc_iter = 2600; jumps = 5;
nu = 6; Psi = diag(1,p); lambda = 0.01; mu = rep(0, p)

cell_res1 <- run_mcmc(df, p,
                      burn_in = burn_in, jumps = jumps, mcmc_iter = mcmc_iter,
                      alpha = alpha, beta = beta, S_init = 1, r_init = 1,
                      nu = nu, Psi = Psi, lambda = lambda, mu = mu, seed = NULL,
                      file_pre = "/mNDP Results/", file_post = "_cell_res1_v2")
saveRDS(cell_res1, here("mNDP Results", "cell_res1_v2.RDS"))
rm(cell_res1)

cell_res2 <- run_mcmc(df, p,
                      burn_in = burn_in, jumps = jumps, mcmc_iter = mcmc_iter,
                      alpha = alpha, beta = beta, S_init = 12, r_init = 2,
                      nu = nu, Psi = Psi, lambda = lambda, mu = mu, seed = NULL,
                      file_pre = "/mNDP Results/", file_post = "_cell_res2_v2")
saveRDS(cell_res2, here("mNDP Results", "cell_res2_v2.RDS"))
rm(cell_res2)

cell_res3 <- run_mcmc(df, p,
                      burn_in = burn_in, jumps = jumps, mcmc_iter = mcmc_iter,
                      alpha = alpha, beta = beta, S_init = 5, r_init = 3,
                      nu = nu, Psi = Psi, lambda = lambda, mu = mu, seed = NULL,
                      file_pre = "/mNDP Results/", file_post = "_cell_res3_v2")
saveRDS(cell_res3, here("mNDP Results", "cell_res3_v2.RDS"))
rm(cell_res3)

####### Tissue Dataset ##########
df <- as.matrix(tissue_df); p <- dim(df)[2] - 1
# remove groups with fewer than 2 observations

# Set parameters as per their paper, adjust nu/psi to wishart as opposed to gamma #
alpha = 0.5; beta = 2; burn_in = 2000; mcmc_iter = 1600; jumps = 5;
nu = 6; Psi = diag(1,p); lambda = 0.01; mu = rep(0, p)

tissue_res1 <- run_mcmc(df, p,
                        burn_in = burn_in, jumps = jumps, mcmc_iter = mcmc_iter,
                        alpha = alpha, beta = beta, S_init = 1, r_init = 1,
                        nu = nu, Psi = Psi, lambda = lambda, mu = mu, seed = NULL,debug = F,
                        file_pre = "/mNDP Results/", file_post = "_tissue_res1_v2")
saveRDS(tissue_res1, here("mNDP Results", "tissue_res1_v2.RDS"))
rm(tissue_res1)

tissue_res2 <- run_mcmc(df, p,
                        burn_in = burn_in, jumps = jumps, mcmc_iter = mcmc_iter,
                        alpha = alpha, beta = beta, S_init = 12, r_init = 2,
                        nu = nu, Psi = Psi, lambda = lambda, mu = mu, seed = NULL,
                        file_pre = "/mNDP Results/", file_post = "_tissue_res2_v2")
saveRDS(tissue_res2, here("mNDP Results", "tissue_res2_v2.RDS"))
rm(tissue_res2)

tissue_res3 <- run_mcmc(df, p,
                        burn_in = burn_in, jumps = jumps, mcmc_iter = mcmc_iter,
                        alpha = alpha, beta = beta, S_init = 5, r_init = 3,
                        nu = nu, Psi = Psi, lambda = lambda, mu = mu, seed = NULL,
                        file_pre = "/mNDP Results/", file_post = "_tissue_res3_v2")
saveRDS(tissue_res3, here("mNDP Results", "tissue_res3_v2.RDS"))
rm(tissue_res3)

# ######### Results Cell df#########
# # Cell
df <- readRDS(here("RNA Splicing Data", "cell group df.RDS"))
res1 <- readRDS(here("mNDP Results", "cell_res1_v2.RDS"))
res2 <- readRDS(here("mNDP Results", "cell_res2_v2.RDS"))
res3 <- readRDS(here("mNDP Results", "cell_res3_v2.RDS"))

# Gelman diagnostic
library(coda)
mcmc.res1 <- mcmc(data= res1$total_l_p, thin = 5)
mcmc.res2 <- mcmc(data= res2$total_l_p, thin = 5)
mcmc.res3 <- mcmc(data= res3$total_l_p, thin = 5)

mcmc_list <- list(mcmc.res1, mcmc.res2, mcmc.res3)

gelman.diag(mcmc_list, confidence = 0.95, transform=TRUE, autoburnin=TRUE,
            multivariate=FALSE)
#
p1 <- plot_results(df, res1)
p2 <- plot_results(df, res2)
p3 <- plot_results(df, res3)
gridExtra::grid.arrange(p1,p2, p3)
#

# Log prob plot
idx <- seq(0,length(res1$total_l_p), by = 5)
ylim <- c(min(c(res1$total_l_p[idx], res2$total_l_p[idx],res3$total_l_p[idx])),
          max(c(res1$total_l_p[idx], res2$total_l_p[idx],res3$total_l_p[idx])))
plot(res2$total_l_p[idx], ylim = ylim, type = "l", col = "red")
lines(res1$total_l_p[idx], col = "blue")
lines(res3$total_l_p[idx], col = "green")
abline(v = burn_in/5)

# Total K - already doesn't include burn_in
ylim <- c(min(c(res1$total_K, res2$total_K,res3$total_K)),
          max(c(res1$total_K, res2$total_K,res3$total_K)))
plot(res2$total_K, ylim = ylim, type = "l", col = "red")
lines(res1$total_K, col = "blue")
lines(res3$total_K, col = "green")

# ######### Results Tissue #########
df <- readRDS(here("RNA Splicing Data", "tissue group df.RDS"))
res1 <- readRDS(here("mNDP Results", "tissue_res1_v2.RDS"))
res2 <- readRDS(here("mNDP Results", "tissue_res2_v2.RDS"))
res3 <- readRDS(here("mNDP Results", "tissue_res3_v2.RDS"))

# # Gelman diagnostic
library(coda)
mcmc.res1 <- mcmc(data= res1$total_l_p, thin = 5)
mcmc.res2 <- mcmc(data= res2$total_l_p, thin = 5)
mcmc.res3 <- mcmc(data= res3$total_l_p, thin = 5)

mcmc_list <- list(mcmc.res1, mcmc.res2, mcmc.res3)

gelman.diag(mcmc_list, confidence = 0.95, transform=TRUE, autoburnin=TRUE,
            multivariate=FALSE)


p1 <- plot_results(df, res1)
p2 <- plot_results(df, res2)
p3 <- plot_results(df, res3)

gridExtra::grid.arrange(p1,p2, p3)


ggplot(df, aes(x = X1, y = X2, color = factor(Cell_sys_num)))



####### modify alpha beta ##########
df <- as.matrix(cell_df); p <- dim(df)[2] - 1

# Set parameters as per their paper, adjust nu/psi to wishart as opposed to gamma #
alpha = 0.5; beta = 0.3; burn_in = 1000; mcmc_iter = 1000; jumps = 5;
nu = 6; Psi = diag(1,p); lambda = 0.01; mu = rep(0, p)

cell_res1 <- run_mcmc(df, p,
                      burn_in = burn_in, jumps = jumps, mcmc_iter = mcmc_iter,
                      alpha = alpha, beta = beta, S_init = 1, r_init = 1,
                      nu = nu, Psi = Psi, lambda = lambda, mu = mu, seed = NULL,
                      file_pre = "/mNDP Results/", file_post = "_cell_res1_v3")
saveRDS(cell_res1, here("mNDP Results", "cell_res1_v3.RDS"))
plot_results(df, cell_res1)
rm(cell_res1)

# Set parameters as per their paper, adjust nu/psi to wishart as opposed to gamma #
alpha = 0.5; beta = 1; burn_in = 100; mcmc_iter = 100; jumps = 5;
nu = 10; Psi = diag(1,p); lambda = 10; mu = rep(0, p)

cell_res1 <- run_mcmc(df, p,
                      burn_in = burn_in, jumps = jumps, mcmc_iter = mcmc_iter,
                      alpha = alpha, beta = beta, S_init = 1, r_init = 1,
                      nu = nu, Psi = Psi, lambda = lambda, mu = mu, seed = NULL,
                      file_pre = "/mNDP Results/", file_post = "_cell_res1_v4")
saveRDS(cell_res1, here("mNDP Results", "cell_res1_v4.RDS"))
plot_results(df, cell_res1)
