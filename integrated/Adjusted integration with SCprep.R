# Load libraries
library(Seurat)
library(monocle3)
library(ggplot2)
library(SeuratWrappers)
library(SummarizedExperiment)
library(magrittr)
library(SeuratData)

set.seed(42)

#read LRP seurat object
LRP.dataset <- readRDS("C:/Lateral root Primordia R project/Visualization tool/LRP.rds")

#read and converting the monocle object to seurat
# Specify the path of your file
file_path <- "GSE158761_all_cells_labeled.rds.gz"

# Load the compressed .rds file
con <- gzcon(gzfile(file_path, "rb"))
monocle_object <- readRDS(con)
close(con)

# Check the class of the loaded data
class(monocle_object)

# Get the counts matrix
counts <- SummarizedExperiment::assay(monocle_object)

# Extract the metadata
meta.data <- as.data.frame(monocle_object@colData)

# Create a Seurat object
mono_seurat <- CreateSeuratObject(counts = counts, meta.data = meta.data)

#check class
class(mono_seurat)

#read txt file to convert to seuratobject
file_path <- "GSE161970_SCO_LRP.txt.gz"

# read the compressed .txt file
new_data <- read.table(file_path, header = TRUE, sep="\t")

# Make the "Locus" column the row names of the data frame
rownames(new_data) <- new_data$Locus

# Remove the "Locus" column from the data frame
new_data$Locus <- NULL

# Create the Seurat object after placing data in new row
new_seurat <- CreateSeuratObject(counts = new_data, project = "GSE161970")

# Create a list of the datasets
datasets <- list(LRP.dataset, mono_seurat, new_seurat)

# SCTransform on each dataset in the list
datasets <- lapply(X = datasets, FUN = SCTransform, verbose = FALSE)

# select integration features
features <- SelectIntegrationFeatures(object.list = datasets)

# prep step for integration
datasets <- lapply(X = datasets, FUN = function(x) {
  x <- PrepSCTIntegration(object = x, anchor.features = features)
  return(x)
})

# downstream integration steps
anchors <- FindIntegrationAnchors(
  object.list = datasets,
  normalization.method = "SCT",
  anchor.features = features
)

integrated_data <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

# Set the default assay to "integrated"
DefaultAssay(integrated_data) <- 'integrated'

# Run PCA
integrated_data <- RunPCA(integrated_data, verbose = FALSE)

# Run UMAP for visualization
integrated_data <- RunUMAP(integrated_data, dims = 1:30, verbose = FALSE)

# Find neighbors
integrated_data <- FindNeighbors(integrated_data, dims = 1:30, verbose = FALSE)

# Identify cell clusters
integrated_data <- FindClusters(integrated_data, resolution = 0.5, verbose = FALSE)

# Save the processed data
saveRDS(integrated_data, 'Integrated.processed2.rds')