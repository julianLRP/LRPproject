# Load gene_info_data
gene_info_data <- read.csv("gene_info.csv", stringsAsFactors = FALSE)

# Define your list of resolutions
resolutions <- c(0.2, 0.5, 1, 1.5, 2)

# Iterate through unique cell types
unique_cell_types <- unique(gene_info_data$cell_type)

# Gather all genes according to cell type order
all_genes <- unlist(lapply(unique_cell_types, function(ct) {
  gene_info_data[gene_info_data$cell_type == ct,]$gene_id
}))

# Filter out genes that do not exist in the dataset
all_genes <- all_genes[all_genes %in% rownames(processed_LRP[["integrated"]]@counts)]

# Remove duplicates
all_genes <- unique(all_genes)

# Get the start index of each cell type
start_indices <- sapply(unique_cell_types, function(ct) {
  match(gene_info_data[gene_info_data$cell_type == ct,]$gene_id[1], all_genes)
})

# Create a data frame for cell type annotations
cell_type_annotations <- data.frame(
  cell_type = unique_cell_types,
  start_index = start_indices
)

# Sort this data frame by start_index
cell_type_annotations <- cell_type_annotations[order(cell_type_annotations$start_index), ]

# Outer loop iterating over resolutions
for(res in resolutions) {
  
  # Find neighbors for the current resolution
  processed_LRP <- FindNeighbors(processed_LRP, dims = 1:50, verbose = FALSE)
  
  # Create clusters for the current resolution
  processed_LRP <- FindClusters(processed_LRP, resolution = res, verbose = FALSE)
  
  # Define directory name
  dir_name <- paste0("Dotplots/Resolution_", res)
  
  # Create directory if it doesn't exist
  if (!dir.exists(dir_name)) {
    dir.create(dir_name, recursive = TRUE)
  }
  
  # Define the PNG file name
  png_file_name <- file.path(dir_name, paste0("Resolution_", res, "_AllGenes.png"))
  
  # Check if the file already exists, and if not, generate the plots and save the file
  if (!file.exists(png_file_name)) {
    # Create a dot plot for the gene
    dot_plot <- DotPlot(processed_LRP, features = all_genes, assay = "integrated") +
      scale_color_gradientn(colors = c('lightgray','yellow','red','darkred')) +
      theme(plot.background = element_rect(fill = "white"),  # Make the background white
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust =1, size = 8)) +  # Rotate X axis gene names 90 degrees
      ggtitle(paste0("Resolution: ", res))
    
    # Add vertical lines or text to denote cell types
    for(i in 1:(nrow(cell_type_annotations) - 1)) {
      dot_plot <- dot_plot + 
        geom_vline(xintercept = cell_type_annotations$start_index[i], color = "black", linetype = "dashed") +
        annotate("text", x = (cell_type_annotations$start_index[i] + cell_type_annotations$start_index[i+1]) / 2, y = 1, label = cell_type_annotations$cell_type[i])
    }
    
    # Add the last cell type
    dot_plot <- dot_plot + 
      annotate("text", x = (cell_type_annotations$start_index[nrow(cell_type_annotations)] + length(all_genes)) / 2, y = 1, label = cell_type_annotations$cell_type[nrow(cell_type_annotations)])
    
    # Save the dot plot as a PNG file, specifying a larger width
    ggsave(filename = png_file_name, plot = dot_plot, width = 40, height = 16)  # Double the width
  }
}
