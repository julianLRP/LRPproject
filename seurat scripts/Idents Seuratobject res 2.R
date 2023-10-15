library(Seurat)
library(slingshot)
library(SingleCellExperiment)
library(ggplot2)

set.seed(42)

# Load your processed LRP dataset
processed_LRP <- readRDS('Integrated.processed2.rds')

# Set the resolution to 2
res <- 2

# Create clusters for the set resolution
processed_LRP <- FindClusters(processed_LRP, resolution = res)

# Rename the clusters
processed_LRP <- RenameIdents(
  object = processed_LRP,
  `0` = "Pericycle Phloem pole and or Procambium",
  `1` = "Pericycle xylem pole and or Phloem",
  `2` = "Procambium and or stele",
  `3` = "Pericycle of stage unknown 1",
  `4` = "Pericycle of stage unknown 2",
  `5` = "Pericycle phloem pole",
  `6` = "Unknown",
  `7` = "Endodermis and or artichoblast",
  `8` = "Phloem Perycicle",
  `9` = "Metaxylem Pericycle",
  `10` = "Columella LRC",
  `11` = "LRC",
  `12` = "Endodermis epidermis Metaxylem initials",
  `13` = "Trichoblast",
  `14` = "Procambium",
  `15` = "Pericycle of stage unknown 3",
  `16` = "Endodermis more mature",
  `17` = "Procambium 2",
  `18` = "Initials 1",
  `19` = "Procambium 3",
  `20` = "Cortex 2",
  `21` = "Casparian strip cluster",
  `22` = "Initials 2",
  `23` = "Outer layer of LRP",
  `24` = "Pericycle of stage unknown 5",
  `25` = "Metaxylem QC",
  `26` = "Unknown 2",
  `27` = "Dividing cells",
  `28` = "Cortex",
  `29` = "Outer layer of LRP",
  `30` = "Trichoblast",
  `31` = "LRC",
  `32` = "Protoxylem",
  `33` = "Unknown 4",
  `34` = "Unknown 5",
  `35` = "Potential initials"
)

# Save the processed LRP dataset
saveRDS(processed_LRP, 'LRP.processed2.rds')

gene_info_data <- read.csv("gene_info.csv", stringsAsFactors = FALSE)

# Set the default assay to RNA (if you have it in your object)
DefaultAssay(processed_LRP) <- "RNA"

# FindNeighbors step before entering the loop over resolutions
processed_LRP <- FindNeighbors(processed_LRP, dims = 1:20, verbose = FALSE)

# Define your list of resolutions
resolutions <- c(2)

# Iterate through unique cell types
unique_cell_types <- unique(gene_info_data$cell_type)

# Gather all genes according to cell type order
all_genes <- unlist(lapply(unique_cell_types, function(ct) {
  gene_info_data[gene_info_data$cell_type == ct,]$gene_id
}))

# Filter out genes that do not exist in the dataset
all_genes <- all_genes[all_genes %in% rownames(processed_LRP[["RNA"]]@counts)]

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
  dot_plot <- DotPlot(processed_LRP, features = all_genes) +
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