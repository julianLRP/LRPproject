install.packages("Matrix")
install.packages("Seurat")
install.packages("dplyr", type = "source")
install.packages(" ggplot2")
install.packages("devtools")
install.packages("dittoSeq")
library(devtools)
install_github("dtm2451/dittoSeq")
install.packages("tidyverse")
install.packages("conflicted")
install.packages("patchwork")

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(dittoSeq)
library(tidyverse)
library(conflicted)
library(patchwork)

set.seed(42)

# Updated path for the LRP RDS file
LRP.dataset <- readRDS("C:/Lateral root Primordia R project/Visualization tool/LRP.rds")

# Process the LRP dataset
processed_LRP <- SCTransform(LRP.dataset)
processed_LRP <- RunPCA(processed_LRP, verbose = FALSE)
processed_LRP <- RunUMAP(processed_LRP, dims = 1:50, verbose = FALSE)
processed_LRP <- FindNeighbors(processed_LRP, dims = 1:50, verbose = FALSE)
processed_LRP <- FindClusters(processed_LRP, resolution = 1, verbose = FALSE)

# Save the processed LRP dataset
saveRDS(processed_LRP, 'LRP.processed.rds')

# Load the processed LRP dataset
processed_LRP <- readRDS("C:/Lateral root Primordia R project/Visualization tool/LRP.processed.rds")

# Calculate the percentage of mitochondrial and chloroplast transcripts for pre-SCT dataset
LRP.dataset[["percent.mt.preSCT"]] <- PercentageFeatureSet(LRP.dataset, pattern = "^ATM")
LRP.dataset[["percent.ct.preSCT"]] <- PercentageFeatureSet(LRP.dataset, pattern = "^ATC")

# Calculate the percentage of mitochondrial and chloroplast transcripts for post-SCT dataset
processed_LRP[["percent.mt.postSCT"]] <- PercentageFeatureSet(processed_LRP, pattern = "^ATM")
processed_LRP[["percent.ct.postSCT"]] <- PercentageFeatureSet(processed_LRP, pattern = "^ATC")

# Create violin plots for pre-SCT dataset
violin_preSCT <- VlnPlot(
  object = LRP.dataset,
  features = c("percent.ct.preSCT", "percent.mt.preSCT"),
  ncol = 2
) + ggtitle("Pre-SCTransform")

# Create violin plots for post-SCT dataset
violin_postSCT <- VlnPlot(
  object = processed_LRP,
  features = c("percent.ct.postSCT", "percent.mt.postSCT"),
  ncol = 2
) + ggtitle("Post-SCTransform")

# Combine and display violin plots for pre-SCT and post-SCT datasets
combined_plot <- violin_preSCT / violin_postSCT
combined_plot

DimPlot(processed_LRP, label = TRUE)

# Identify marker genes specific to each cluster
all_markers <- FindAllMarkers(processed_LRP, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

top_markers <- all_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

print(top_markers)
top_marker_genes <- unique(top_markers$gene)

print(top_markers, n = 210)

# Make large heatmap
pdf("heatmap_large.pdf", width = 20, height = 25)  
DoHeatmap(processed_LRP, features = top_marker_genes) + NoLegend()
