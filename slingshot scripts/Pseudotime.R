  library(monocle3)
  library(magrittr)
  library(Seurat)
  library(SeuratData)
  library(SeuratWrappers)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(dittoSeq)
  library(tidyverse)
  library(conflicted)
  library(patchwork)
  library(gridExtra)
  library(cowplot)
  library(stringr)
  library(slingshot)
  library(SingleCellExperiment)

  
#read processed dataset
  processed_LRP <- readRDS("C:/Lateral root Primordia R project/Visualization tool/LRP.processed.rds")
  
  # Extract cell clusters and UMAP coordinates from Seurat object
  umap_coords <- processed_LRP[["umap"]]@cell.embeddings
  cl <- processed_LRP$RNA_snn_res.1  # cluster identities
  
  # Prepare a SingleCellExperiment object
  sce <- SingleCellExperiment(assays = list(counts = t(umap_coords)))
  
  # Make sure cluster identities match cells in sce
  cl_aligned <- factor(cl[colnames(sce)])
  
  # Assign cluster identities to colData of sce
  colData(sce)$Cluster <- cl_aligned
  
  # Add UMAP coordinates as a reduced dimension
  reducedDims(sce)$UMAP <- umap_coords
  
  # Run slingshot
  sce <- slingshot(sce, clusterLabels = 'Cluster', reducedDim = 'UMAP')
  
  # Convert factor to numeric for color assignment
  cluster_numeric <- as.numeric(sce$Cluster)
  
  # Identify the unique cluster names
  unique_clusters <- unique(colData(sce)$Cluster)
  
  # Plotting the lineages
  plot(reducedDims(sce)$UMAP, col = rainbow(max(cluster_numeric))[cluster_numeric],
       main = 'Slingshot Lineages')
  lines(SlingshotDataSet(sce), lwd = 2, type = 'lineages')
  
# Add cluster names to the plot
for (cluster in unique_clusters) {
  cluster_coords <- reducedDims(sce)$UMAP[colData(sce)$Cluster == cluster, ]
  cluster_center <- apply(cluster_coords, 2, mean)
  text(cluster_center[1], cluster_center[2], labels = cluster, cex = 1.2, pos = 1)
}

  # Save the modified 'sce' object with pseudotime as an RDS file
  saveRDS(sce, "LRP_with_pseudotime.rds")
