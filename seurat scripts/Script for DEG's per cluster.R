library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(dplyr)
library(tidyr)
library(ggplot2)

# Initialize an empty list to store the results
all.markers <- list()

# Loop over the cluster names
for(cluster in unique(Idents(processed_LRP))) {
  # Find markers for the current cluster
  markers <- FindMarkers(processed_LRP, ident.1 = cluster)
  
  # Store the results in the list
  all.markers[[cluster]] <- markers
}

# Loop over the list
for(cluster in names(all.markers)) {
  # Replace slashes in the cluster name with underscores
  safe_cluster_name <- gsub("/", "_", cluster)
  
  # Write the data frame to a CSV file
  write.csv(all.markers[[cluster]], file = paste0(safe_cluster_name, "_markers.csv"))
}