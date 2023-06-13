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

# Pre-processing steps
processed_LRP <- SCTransform(LRP.dataset)
DefaultAssay(processed_LRP) <- 'RNA'
processed_LRP <- RunPCA(processed_LRP, verbose = FALSE)
processed_LRP <- RunUMAP(processed_LRP, dims = 1:50, verbose = FALSE)

# FindNeighbors step before entering the loop over resolutions
processed_LRP <- FindNeighbors(processed_LRP, dims = 1:50, verbose = FALSE)

# Define your list of resolutions
resolutions <- c(0.2, 0.5, 1, 1.5, 2)

# Outer loop iterating over resolutions
for(res in resolutions) {
  # Create clusters for the current resolution
  processed_LRP <- FindClusters(processed_LRP, resolution = res)
  
  # Iterate through the gene_info_data data frame and create plots for each gene
  for (i in 1:nrow(gene_info_data)) {
    cell_type <- gene_info_data$cell_type[i]
    gene_id <- gene_info_data$gene_id[i]
    gene_name <- gene_info_data$gene_name[i]
    gene_description <- get_gene_description(gene_id)
    
    # Check if the gene exists in the dataset
    if(!gene_id %in% rownames(processed_LRP[["RNA"]]@counts)) {
      print(paste("Skipping gene_id", gene_id, "as it's not present in the dataset."))  # Print message to the console
      next  # Skip to the next iteration of the loop
    }
    
    # Define directory name
    dir_name <- paste0("Resolution_", res, "/", cell_type)
    
    # Create directory if it doesn't exist
    if (!dir.exists(dir_name)) {
      dir.create(dir_name, recursive = TRUE)
    }
    
    # Define the PNG file name
    png_file_name <- file.path(dir_name, paste0("Resolution_", res, "_", cell_type, "_genes_", gene_id, "_", gene_name, ".png"))
    
    # Check if the file already exists, and if not, generate the plots and save the file
    if (!file.exists(png_file_name)) {
      # Create a feature plot for the gene
      feature_plot <- FeaturePlot(processed_LRP, features = gene_id, order = T, label = T, pt.size = 3) +
        scale_color_gradientn(colors = c('lightgray','yellow','red','darkred')) +
        ggtitle(paste0(gene_id, " - ", gene_name))
      
      # Add the gene of interest to the metadata
      processed_LRP[["gene_of_interest"]] <- ifelse(processed_LRP[["RNA"]]@counts[gene_id,] > 0, gene_id, "Other")
      
      # Create a dim plot for the gene
      color_vector <- setNames(c(rgb(1,0,0,alpha=0.5), rgb(0.5,0.5,0.5,alpha=0.5)), c(gene_id, "Other"))
      dim_plot <- DimPlot(processed_LRP, group.by = "gene_of_interest", pt.size = 3, cols = color_vector) +
        ggtitle(paste0(gene_id, " - ", gene_name))
      
      # Create a description text box
      wrapped_gene_description <- str_wrap(gene_description, width = 80)  # Adjust the width as needed
      description_lines <- str_split(wrapped_gene_description, "\n", simplify = FALSE)[[1]]
      line_count <- length(description_lines)
      
      # Create a data frame with the lines and their corresponding y-positions
      y_increment <- 1 / (line_count + 1)  # Adjust this value to increase or decrease the space between lines
      description_df <- data.frame(
        x = rep(0, line_count),
        y = seq(from = 1 - y_increment, to = 0, length.out = line_count),
        text = description_lines
      )
      
      # Calculate font size based on the number of lines
      if (line_count <= 10) {
        font_size <- 5
      } else if (line_count <= 20) {
        font_size <- 4
      } else {
        font_size <- 3
      }
      
      description_text <- ggplot(data = description_df) +
        geom_text(aes(x = x, y = y, label = text), hjust = 0, vjust = 0, size = font_size, family = "mono") +  # Use the calculated font_size
        theme_void() +
        theme(plot.margin = margin(t = 5, r = 5, b = 5, l = 5, unit = "pt"),
              text = element_text(family = "mono")) +
        coord_cartesian(clip = 'off', xlim = c(0, 1), ylim = c(0, 1))
      
      # Combine the plots and the description text box
      combined_plot <- (feature_plot / dim_plot / description_text) +
        patchwork::plot_layout(ncol = 1, heights = c(2, 2, 1))
      
      # Save the combined plot as a PNG file, specifying a larger width and height
      ggsave(filename = png_file_name, plot = combined_plot, width = 12, height = 16)
      
      # Remove the gene of interest from the metadata to prevent confusion in the next iteration
      processed_LRP[["gene_of_interest"]] <- NULL
    }
  }
}