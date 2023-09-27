library(Seurat)
library(ggplot2)

# Load your processed LRP dataset
processed_LRP <- readRDS('LRP.processed2.rds')

plot_cluster <- function(seurat_obj, cluster_id) {
  # Fetch the current identities from the Seurat object
  current_idents <- Idents(seurat_obj)
  
  # Create a custom vector of colors: Set your chosen cluster color, and the rest to gray
  custom_colors <- rep("lightgray", length(unique(current_idents)))
  names(custom_colors) <- unique(current_idents)
  custom_colors[as.character(cluster_id)] <- "red"
  
  # Use FetchData to get the UMAP coordinates
  umap_data <- FetchData(seurat_obj, vars = c("UMAP_1", "UMAP_2"))
  umap_data$cluster <- current_idents
  
  plot <- ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2, color = cluster)) +
    geom_point() +
    scale_color_manual(values = custom_colors) +
    labs(title = paste("Cluster", cluster_id)) +
    theme_minimal()
  
  return(plot)
}

# Loop through each cluster and plot:
unique_clusters <- unique(Idents(processed_LRP))
for (cluster in unique_clusters) {
  print(plot_cluster(processed_LRP, cluster))
}


#or use this to save the images
plot_cluster <- function(seurat_obj, cluster_id) {
  # Fetch the current identities from the Seurat object
  current_idents <- Idents(seurat_obj)
  
  # Create a custom vector of colors: Set your chosen cluster color, and the rest to gray
  custom_colors <- rep("lightgray", length(unique(current_idents)))
  names(custom_colors) <- unique(current_idents)
  custom_colors[as.character(cluster_id)] <- "red"
  
  # Use FetchData to get the UMAP coordinates
  umap_data <- FetchData(seurat_obj, vars = c("UMAP_1", "UMAP_2"))
  umap_data$cluster <- current_idents
  
  plot <- ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2, color = cluster)) +
    geom_point() +
    scale_color_manual(values = custom_colors) +
    labs(title = paste("Cluster", cluster_id)) +
    theme_minimal() +
    theme(legend.position = "none")  # Omit the legend
  
  # Save the plot as an image
  file_name <- paste0("Cluster_", cluster_id, ".png")
  ggsave(filename = file_name, plot = plot, width = 6, height = 4)
  
  return(plot)
}

# Loop through each cluster and plot (and save):
unique_clusters <- unique(Idents(processed_LRP))
for (cluster in unique_clusters) {
  plot_cluster(processed_LRP, cluster)
}