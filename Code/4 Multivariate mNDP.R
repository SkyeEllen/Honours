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


####### Run All ##########
cell_df <- as.matrix(cell_df); cell_p <- dim(cell_df)[2] - 1
tissue_df <- as.matrix(tissue_df); tissue_p <- dim(tissue_df)[2]-1
# Set parameters as per their paper, adjust nu/psi to wishart as opposed to gamma #
alpha = 0.5; beta = 1; burn_in = 2000; mcmc_iter = 1000; jumps = 5;

S_init <- c(1, 3, 5, 7, 9); r_init <- c(1, 4, 2, 3, 1)

for(i in 1:5){
  nu = 6; Psi = diag(1,cell_p); lambda = 0.01; mu = rep(0, cell_p)
  res <- run_mcmc(cell_df, cell_p,
                  burn_in = burn_in, jumps = jumps, mcmc_iter = mcmc_iter,
                  alpha = alpha, beta = beta, S_init = S_init[i], r_init = r_init[i],
                  nu = nu, Psi = Psi, lambda = lambda, mu = mu, seed = NULL,
                  file_pre = "/mNDP New Results 2/", file_post = paste("_cell_res", i, sep = ""))
  saveRDS(res, here("mNDP New Results 2", paste("cell_res", i, ".RDS", sep = "")))

  nu = 6; Psi = diag(1,tissue_p); lambda = 0.01; mu = rep(0, tissue_p)
  res <- run_mcmc(tissue_df, tissue_p,
                  burn_in = burn_in, jumps = jumps, mcmc_iter = mcmc_iter,
                  alpha = alpha, beta = beta, S_init = S_init[i], r_init = r_init[i],
                  nu = nu, Psi = Psi, lambda = lambda, mu = mu, seed = NULL,
                  file_pre = "/mNDP New Results 2/", file_post = paste("_tissue_res", i, sep = ""))
  saveRDS(res, here("mNDP New Results 2", paste("tissue_res", i, ".RDS", sep = "")))
}


# ######### Results Cell df#########
# # Cell
df <- readRDS(here("RNA Splicing Data", "cell group df.RDS"))
res <- list(); for(i in 1:max_res){res[[i]] <- readRDS(here("mNDP New Results", paste("cell_res", i, ".RDS", sep = ""))) }
mcmc_list <- list(); idx <- seq(burn_in + 1, length(res[[1]]$total_l_p))
for(i in 1:max_res) {mcmc_list[[i]] <- mcmc(res[[i]]$total_l_p[idx])}
gelman.diag(mcmc_list, autoburnin = F)

plots <- list(); for(i in 1:max_res) plots[[i]] <- plot_results(df, res[[i]])
gridExtra::grid.arrange(grobs = plots)
# All only have one observational cluster

# Log prob plot
idx <- seq(0,length(res[[i]]$total_l_p), by = 5)
min <- 0; for(i in 1:max_res) min <- min(min, res[[i]]$total_l_p)
max <- min; for(i in 1:max_res) max <- max(max, res[[i]]$total_l_p)
cols <- RColorBrewer::brewer.pal(5, "Set2")
plot(res[[1]]$total_l_p[idx], ylim = c(min, max), type = "l", col = cols[1])
for(i in 2:max_res) lines(res[[i]]$total_l_p[idx], col = cols[i])
abline(v = burn_in/5, lty = 3)

# Total K
min <- 0; for(i in 1:max_res) min <- min(min, res[[i]]$total_K)
max <- min; for(i in 1:max_res) max <- max(max, res[[i]]$total_K)
cols <- RColorBrewer::brewer.pal(5, "Set2")
plot(res[[1]]$total_K[idx], ylim = c(min, max), type = "l", col = cols[1])
for(i in 2:max_res) lines(res[[i]]$total_K[idx], col = cols[i])
abline(v = burn_in/5, lty = 3)


# ######### Results Tissue #########
# # Cell
df <- readRDS(here("RNA Splicing Data", "tissue group df.RDS"))
res <- list(); for(i in 1:max_res){res[[i]] <- readRDS(here("mNDP New Results", paste("tissue_res", i, ".RDS", sep = ""))) }
mcmc_list <- list(); idx <- seq(burn_in + 1, length(res[[1]]$total_l_p))
for(i in 1:max_res) {mcmc_list[[i]] <- mcmc(res[[i]]$total_l_p[idx])}
gelman.diag(mcmc_list, autoburnin = F)

plots <- list(); for(i in 1:max_res) plots[[i]] <- plot_results(df, res[[i]])
gridExtra::grid.arrange(grobs = plots)
# All only have one observational cluster

# Log prob plot
idx <- seq(0,length(res[[i]]$total_l_p), by = 5)
min <- 0; for(i in 1:max_res) min <- min(min, res[[i]]$total_l_p)
max <- min; for(i in 1:max_res) max <- max(max, res[[i]]$total_l_p)
cols <- RColorBrewer::brewer.pal(5, "Set2")
plot(res[[1]]$total_l_p[idx], ylim = c(min, max), type = "l", col = cols[1])
for(i in 2:max_res) lines(res[[i]]$total_l_p[idx], col = cols[i])
abline(v = burn_in/5, lty = 3)

# Total K
min <- 0; for(i in 1:max_res) min <- min(min, res[[i]]$total_K)
max <- min; for(i in 1:max_res) max <- max(max, res[[i]]$total_K)
cols <- RColorBrewer::brewer.pal(5, "Set2")
plot(res[[1]]$total_K[idx], ylim = c(min, max), type = "l", col = cols[1])
for(i in 2:max_res) lines(res[[i]]$total_K[idx], col = cols[i])
abline(v = burn_in/5, lty = 3)



