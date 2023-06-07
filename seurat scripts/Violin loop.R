for (i in 1:nrow(gene_info_data)) {
  cell_type <- gene_info_data$cell_type[i]
  gene_id <- gene_info_data$gene_id[i]
  gene_name <- gene_info_data$gene_name[i]
  gene_description <- get_gene_description(gene_id)
  
  # Define directory name
  dir_name <- cell_type
  
  # Create directory if it doesn't exist
  if (!dir.exists(dir_name)) {
    dir.create(dir_name)
  }
  
  # Define the PNG file name for the violin plot
  png_file_name_violin <- file.path(dir_name, paste0(cell_type, "_violin_genes_", gene_id, "_", gene_name, ".png"))
  
  # Check if the file already exists, if not, generate the plots and save the file
  if (!file.exists(png_file_name_violin)) {
    # Create a violin plot for the gene
    violin_plot <- VlnPlot(processed_LRP, features = gene_id, group.by = "seurat_clusters")
    
    # Save the violin plot as a PNG file
    ggsave(filename = png_file_name_violin, plot = violin_plot, width = 12, height = 8)
  }
}