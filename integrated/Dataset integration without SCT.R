# Load libraries
library(Seurat)
library(monocle3)
library(ggplot2)
library(SeuratWrappers)
library(SummarizedExperiment)
library(magrittr)

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

# Calculate percent of mitochondrial and chloroplast genes
LRP.dataset[["percent.mt"]] <- PercentageFeatureSet(LRP.dataset, pattern = "^ATM")
LRP.dataset[["percent.ct"]] <- PercentageFeatureSet(LRP.dataset, pattern = "^ATC")

# Create violin plots
VlnPlot(object = LRP.dataset, features = c("nFeature_RNA", "nCount_RNA","percent.ct","percent.mt"), ncol = 5)

# Filter cells based on QC metrics
LRP.dataset <- subset(LRP.dataset, subset = percent.mt <= 20 & percent.ct <= 20 & nCount_RNA >=500)

# Calculate percent of mitochondrial and chloroplast genes
mono_seurat[["percent.mt"]] <- PercentageFeatureSet(mono_seurat, pattern = "^ATM")
mono_seurat[["percent.ct"]] <- PercentageFeatureSet(mono_seurat, pattern = "^ATC")

# Create violin plots
VlnPlot(object = mono_seurat, features = c("nFeature_RNA", "nCount_RNA","percent.ct","percent.mt"), ncol = 5)

# Filter cells based on QC metrics
mono_seurat <- subset(mono_seurat, subset = percent.mt <= 20 & percent.ct <= 20 & nCount_RNA >=500)

# Calculate percent of mitochondrial and chloroplast genes
new_seurat[["percent.mt"]] <- PercentageFeatureSet(new_seurat, pattern = "^ATM")
new_seurat[["percent.ct"]] <- PercentageFeatureSet(new_seurat, pattern = "^ATC")

# Create violin plots
VlnPlot(object = new_seurat, features = c("nFeature_RNA", "nCount_RNA","percent.ct","percent.mt"), ncol = 5)

# Filter cells based on QC metrics
new_seurat <- subset(new_seurat, subset = percent.mt <= 20 & percent.ct <= 20 & nCount_RNA >=500)

# Normalization
LRP.dataset <- NormalizeData(LRP.dataset)
mono_seurat <- NormalizeData(mono_seurat)
new_seurat <- NormalizeData(new_seurat)

# find integration anchors
anchors <- FindIntegrationAnchors(object.list = list(LRP.dataset, mono_seurat, new_seurat), 
                                  dims = 1:20)

# integrate data
integrated_data <- IntegrateData(anchorset = anchors, dims = 1:20)

# Scale data
integrated_data <- ScaleData(integrated_data, verbose = FALSE)

# Run PCA
integrated_data <- RunPCA(integrated_data, verbose = FALSE)

# Run UMAP for visualization
integrated_data <- RunUMAP(integrated_data, dims = 1:20, verbose = FALSE)

# Find neighbors
integrated_data <- FindNeighbors(integrated_data, dims = 1:20, verbose = FALSE)

# Identify cell clusters
integrated_data <- FindClusters(integrated_data, resolution = 2, verbose = FALSE)

# Save the processed data
saveRDS(integrated_data, 'Integrated.processed2.rds')