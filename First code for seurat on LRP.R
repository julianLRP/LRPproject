#Step 0. Preparatory. Loading libraries.
#Use _packages.R to install all the necessary packages prior loading the libraries

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(dittoSeq)

LRP <- readRDS("C:/Users/Julian/Desktop/GRaduation internship/LRPproject/LRP.rds")

LRP.dataset<- CreateSeuratObject(counts = LRP.data, project = "LRP")
LRP <- SCTransform(LRP)
LRP <- RunPCA(LRP, verbose = FALSE)
LRP <- RunUMAP(LRP, dims = 1:50, verbose = FALSE)
LRP <- FindNeighbors(LRP, dims = 1:50, verbose = FALSE)
LRP <- FindClusters(LRP, resolution = 1, verbose = FALSE)

saveRDS(LRP,'LRP.dataset.rds')

DimPlot(LRP, label = TRUE)

# Identify marker genes specific to each cluster
all_markers <- FindAllMarkers(LRP.dataset, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

top_markers <- all_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

print(top_markers)
top_marker_genes <- unique(top_markers$gene)

print(top_markers,n=210)

#make large heatmap.
pdf("heatmap_large.pdf", width = 20, height = 25)  
DoHeatmap(LRP, features = top_marker_genes) + NoLegend()
