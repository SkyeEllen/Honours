library(umap)       # Calculate UMAP embeddings
library(ggplot2)    # Visualisation of UMAP embeddings / heatmaps
library(progress)   # Progress bar for creating umap plots
library(GGally)     # Facet plot
library(dendextend) # Changing dendrogram aesthetics


################################################################################
# Functions
################################################################################
# Compare various values of n_neighbours and min)dist
create_umap_plots <- function(data, color, seed, n_components){

  # Compare choices of n and d
  n_opts <- seq(from = 5, to = 40,by = 5)
  min_dists <- c(0.001, 0.01,0.05,0.1)

  # Setup plots
  plots <- list(); plots2 <- list()
  i = 1

  # Setup progress bar
  pb <- progress_bar$new(format = "[:bar] :percent ETR: :eta",
                         total = length(n_opts)*length(min_dists),
                         complete = "-", incomplete = " ", current = ">",
                         clear = FALSE, width = 100)

  for(n in n_opts){
    for(d in min_dists){

      # Set UMAP configuration
      umap.config <- umap.defaults
      umap.config$n_neighbors <- n
      umap.config$min_dist <- d
      umap.config$n_components <- n_components

      # Compute UMAP
      set.seed(seed)
      u <- umap(data, config = umap.config, random_state = seed, transform_seed = seed)

      plot.df <- data.frame(u$layout, "col" = color)
      colnames(plot.df)

      # Create plot of results
      p1 <-  ggplot(plot.df) +
        geom_point(aes(x = X1, y = X2, col = factor(col)), show.legend = FALSE) +
        ggtitle(paste("n: ", n, ", d: ", d))
      p2 <- ggplot(plot.df) +
        geom_point(aes(x = X3, y = X1, col = factor(col)), show.legend = FALSE) +
        ggtitle(paste("n: ", n, ", d: ", d))

      # Add to list of plots
      plots[[i]] <- p1
      plots2[[i]] <- p2
      i = i + 1

      pb$tick()
    }
  }

  return(list(plots, plots2))
}

# Create inside a function for reproducibility (when compared to plots from create_umap_plots)
create_umap_single <- function(data, n_neighbors, min_dist, n_components, seed){
  # Set up configuration
  umap.config <- umap.defaults; umap.config$n_neighbors <- n_neighbors
  umap.config$min_dist <- min_dist; umap.config$n_components <- n_components

  # Run umap
  set.seed(seed)
  u <- umap(data, config = umap.config, random_state = seed, transform_seed = seed)
  return(u)
}

create_heatmap <- function(clusters, source){
  hmmat <- table(source, clusters)
  hmmat <- ifelse(hmmat == 0, NA, hmmat)
  h <- heatmap(hmmat, col = hcl.colors((max(hmmat) + 1), "Reds 2", rev = TRUE),
               scale = "row", na.rm = F)
  legend(x="topleft", legend=c(1:max(hmmat)),
         fill=hcl.colors(max(hmmat), "Reds 2", rev = TRUE))
  return(h)
}

create_ggplot_heatmap <- function(clusters, source, xlab = "", ylab = ""){
  hm_df <- expand.grid(cluster = unique(clusters), source = sort(unique(source)))
  count <- as.numeric(table(clusters, source))
  hm_df$count <- ifelse(count != 0, count, NA)
  ggplot(hm_df, aes(x=factor(cluster), y=factor(source), fill=count)) +
    geom_tile(color="white", size = 0.25) +
    scale_fill_gradient(low=RColorBrewer::brewer.pal(3, "Reds")[1],
                        high=RColorBrewer::brewer.pal(3, "Reds")[2],
                        na.value="white") +
    theme_minimal() +
    geom_text(label = ifelse(count > 0, count, "")) +
    labs(x = ifelse(xlab == "", "Cluster", xlab),
         y = ifelse(ylab == "", "System", ylab)) +
    theme(legend.position = "none")
}
