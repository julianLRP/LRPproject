##these scripts probably don't work since it was beyond my technical expertise.

library(Seurat)
library(SingleCellExperiment)
library(slingshot)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(ggrepel)

#perform after Idents Seuratobject res 2
  # Load your processed LRP dataset
  processed_LRP <- readRDS('LRP.processed.rds')
  
  # Define your start cluster
  early_init_cells <- WhichCells(processed_LRP, idents = "Pericycle")
  Idents(processed_LRP, cells = early_init_cells) <- "start.clus"

# Convert your Seurat object to a SingleCellExperiment (SCE) object
sce <- SingleCellExperiment(assays = List(counts = processed_LRP@assays$SCT@counts))
assays(sce)$SCT <- processed_LRP@assays$SCT@data

# Transfer the identities or cluster number
processed_LRP$slingIdent <- Idents(processed_LRP)


# Transfer the embeddings
reducedDims(sce) <- SimpleList(PCA=Embeddings(processed_LRP,reduction = 'pca'), 
                               UMAP = Embeddings(processed_LRP,reduction = 'umap'))
colData(sce)$slingIdent <- processed_LRP$slingIdent

# Run slingshot with PCA or UMAP
sce <- slingshot(sce, reducedDim = 'PCA',clusterLabels = 'slingIdent',start.clus = "start.clus")

# Generate colors for pseudotime
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sce$slingPseudotime_1, breaks=100)]

# Plot PCA colored by pseudotime
plot(reducedDims(sce)$PCA, col = plotcol, pch=16, xlab = "PC1", ylab = "PC2")

# Plot UMAP colored by pseudotime
plot(reducedDims(sce)$UMAP, col = plotcol, pch=16, xlab = "UMAP1", ylab = "UMAP2")

umap_df <- data.frame(UMAP1 = reducedDims(sce)$UMAP[,1], 
                      UMAP2 = reducedDims(sce)$UMAP[,2], 
                      Pseudotime = sce$slingPseudotime_1, 
                      Cluster = processed_LRP@meta.data$slingIdent)

# Create the plot
p <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Pseudotime)) +
  geom_point(pch=16) +
  scale_color_gradientn(colors = colors)

# Calculate the centroids
centroids <- umap_df %>%
  group_by(Cluster) %>%
  summarize(UMAP1 = mean(UMAP1), UMAP2 = mean(UMAP2))

# Add labels with geom_text_repel
p <- p + geom_text_repel(data = centroids, aes(x = UMAP1, y = UMAP2, label = Cluster), size = 3, box.padding = unit(0.35, "lines"), inherit.aes = FALSE)

# Finalize the plot
p <- p + theme_minimal() + labs(color = 'Pseudotime') + 
  theme(legend.position = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
print(p)


