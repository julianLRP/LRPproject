##additional quality control via Elbowplots
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)
library(ggrepel)
library(dplyr)

#Load RDS
integrated_data <- readRDS("Integrated.processed2.rds")

# Set the default assay to RNA (if you have it in your object)
DefaultAssay(integrated_data) <- "RNA"

# Calculate the percentage of mitochondrial and chloroplast transcripts for the integrated dataset
integrated_data[["percent.mt"]] <- PercentageFeatureSet(integrated_data, pattern = "^ATM")
integrated_data[["percent.ct"]] <- PercentageFeatureSet(integrated_data, pattern = "^ATC")

# Create violin plots
violin_plot <- VlnPlot(
  object = integrated_data, 
  features = c("percent.ct", "percent.mt"))

# Save the violin plot
ggsave(
  filename = "combined_violin_plot.png", 
  plot = violin_plot, 
  width = 20, 
  height = 8, 
  dpi = 300
)

# Display the violin plot in the R console
print(violin_plot)


# Display the dimensionality reduction plot with labels
DimPlot(integerated_data , label = TRUE)

if (!dir.exists("Elbowplots")) {
  dir.create("Elbowplots")
}

# Create the Elbow plot object
Elbow_plot <- ElbowPlot(integrated_data)

# Save the plot as a PNG image in the "Elbowplots" directory
ggsave(filename = "Elbowplots/elbow_plot.png", plot = Elbow_plot, dpi = 300)

# Convert to character
integrated_data$orig.ident <- as.character(integrated_data$orig.ident)

# Make the changes
integrated_data$orig.ident[integrated_data$orig.ident == "SeuratProject"] <- "Gala_et_al"
integrated_data$orig.ident[integrated_data$orig.ident == "GSE161970"] <- "Serrano_et_al"

# Convert back to factor
integrated_data$orig.ident <- as.factor(integrated_data$orig.ident)


# Save the object with the renamed idents
saveRDS(integrated_data, file = "Integrated.processed2.rds")

# Assume you have these UMAP coordinates in your Seurat object
umap_df <- data.frame(UMAP1 = integrated_data@reductions$umap@cell.embeddings[,1], 
                      UMAP2 = integrated_data@reductions$umap@cell.embeddings[,2], 
                      orig.ident = integrated_data@meta.data$orig.ident,
                      Cluster = integrated_data@meta.data$seurat_clusters)

# Generate the plot
p <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = orig.ident)) +
  geom_point(pch = 16) +
  scale_color_manual(values = rainbow(length(unique(umap_df$orig.ident))))

# Calculate the centroids
centroids <- umap_df %>%
  group_by(Cluster) %>%
  summarize(UMAP1 = mean(UMAP1), UMAP2 = mean(UMAP2), .groups = 'drop')

# Add labels with geom_text_repel
p <- p + ggrepel::geom_text_repel(data = centroids, aes(x = UMAP1, y = UMAP2, label = Cluster), size = 5, box.padding = unit(0.35, "lines"), inherit.aes = FALSE)

# Finalize the plot
p <- p + theme_minimal() + labs(color = 'orig.ident') + 
  theme(legend.position = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
print(p)
