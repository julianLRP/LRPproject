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
install.packages('BiocManager')
BiocManager::install('limma')
install.packages("gridExtra")
install.packages("SeuratWrappers")
install.packages("magrittr")
install.packages("devtools")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("monocle")
install.packages("devtools", type = "win.binary")
devtools::install_github('satijalab/seurat-data')



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


set.seed(42)

# Updated path for the LRP RDS file
LRP.dataset <- readRDS("C:/Lateral root Primordia R project/Visualization tool/LRP.rds")


# Process the LRP dataset
processed_LRP <- SCTransform(LRP.dataset)
DefaultAssay(processed_LRP)<- 'RNA'
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


# Create violin plots for the pre-SCT dataset
violin_preSCT <- VlnPlot(
  object = LRP.dataset, 
  features = c("percent.ct.preSCT", "percent.mt.preSCT"))
  

violin_postSCT <- VlnPlot(
  object = processed_LRP, 
  features = c("percent.ct.postSCT", "percent.mt.postSCT"))
  

# Combine and display violin plots for pre-SCT and post-SCT datasets
combined_plot_MTandCT <- violin_preSCT / violin_postSCT
ggsave(
  filename = "combined_violin_plot.png", 
  plot = combined_plot_MTandCT, 
  width = 20, 
  height = 8, 
  dpi = 300
)

# Display the combined plot in the R console
print(violin_postSCT)
print(violin_preSCT)
print(combined_plot_MTandCT)


DimPlot(processed_LRP, label = TRUE)

if (!dir.exists("Elbowplots")) {
  dir.create("Elbowplots")
}


# Create the Elbow plot object
Elbow_plot <- ElbowPlot(processed_LRP)

# Save the plot as a PNG image in the "Elbowplots" directory
ggsave(filename = "Elbowplots/elbow_plot.png", plot = Elbow_plot, dpi = 300)

# Identify marker genes specific to each cluster
all_markers <- FindAllMarkers(processed_LRP, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

top_markers <- all_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

print(top_markers)
top_marker_genes <- unique(top_markers$gene)

print(top_markers, n = 210)


# Define gene list
Atrichoblast_gene_list <- list(
  c("AT4G00730", "ANL2"),
  c("AT5G66800", "AT5G66800"),
  c("AT1G68560", "XYL1"),
  c("AT1G31950", "TPS29"),
  c("AT4G25250", "AT4G25250"),
  c("AT1G79840", "GL2"),
  c("AT5G19750", "AT5G19750"),
  c("AT4G00360", "CYP86A2"),
  c("AT4G36250", "ALDH3F1"),
  c("AT5G43040", "AT5G43040"),
  c("AT3G43720", "LTPG2"),
  c("AT1G65310", "XTH17"),
  c("AT2G37260", "WRKY44"),
  c("AT1G68470", "AT1G68470")
)

# Define directory name
dir_name <- "Atrichoblast"

# Create directory if it doesn't exist
if(!dir.exists(dir_name)){
  dir.create(dir_name)
}

for (gene_info in Atrichoblast_gene_list) {
  
  # create a feature plot for the gene
  feature_plot <- FeaturePlot(processed_LRP, features = gene_info[1], order = T, label = T, pt.size = 3) +
    scale_color_gradientn(colors = c('lightgray','yellow','red','darkred')) +
    ggtitle(paste0(gene_info[1], " - ", gene_info[2]))
  
  # Add the gene of interest to the metadata
  processed_LRP[["gene_of_interest"]] <- ifelse(processed_LRP[["RNA"]]@counts[gene_info[1],] > 0, gene_info[1], "Other")
  
  # create a dim plot for the gene
  dim_plot <- DimPlot(processed_LRP, group.by = "gene_of_interest") +
    ggtitle(paste0(gene_info[1], " - ", gene_info[2]))
  
  # Combine the plots into a single grob
  combined_plot <- gridExtra::arrangeGrob(feature_plot, dim_plot, ncol=2)
  
  # Specify the name of the PNG file
  png_file_name <- file.path(dir_name, paste0("Atrichoblast_genes_", gene_info[1], "_", gene_info[2], ".png"))
  
  # Save the combined plot as a PNG file, specifying a larger width and height
  ggsave(filename = png_file_name, plot = combined_plot, width = 15, height = 5)
  
  # Remove the gene of interest from the metadata to prevent confusion in the next iteration
  processed_LRP[["gene_of_interest"]] <- NULL
}

saveRDS(processed_LRP, 'LRP.processed.rds')

#Create a new list of gene IDs and names
# Define gene list
Trichoblast_gene_list <- list(
  c("AT1G03550", "SCAMP2"),
  c("AT2G46380", "AT2G46380"),
  c("AT3G51410", "AT3G51410"),
  c("AT1G33090", "AT1G33090"),
  c("AT3G09330", "AT3G09330"),
  c("AT5G49270", "COBL9"),
  c("AT1G48930", "GH9C1"),
  c("AT1G69240", "MES15"),
  c("AT1G45545", "WEL2"),
  c("AT1G55440", "AT1G55440"),
  c("AT2G39690", "AT2G39690"),
  c("AT1G07795", "AT1G07795"),
  c("AT4G04760", "AT4G04760"),
  c("AT3G60630", "SCL22"),
  c("AT4G04930", "DES-1-LIKE"),
  c("AT1G13950", "ELF5A-1"),
  c("AT2G03220", "FUT1"),
  c("AT4G07960", "CSLC12"),
  c("AT4G25160", "PUB35"),
  c("AT5G07450", "CYCU4-2"),
  c("AT5G65160", "TPR14"),
  c("AT4G09500", "UGT79B7"),
  c("AT3G46270", "AT3G46270"),
  c("AT4G16920", "AT4G16920"),
  c("AT1G07795", "AT1G07795")
)

# Define directory name
dir_name <- "Trichoblast"

# Create directory if it doesn't exist
if(!dir.exists(dir_name)){
  dir.create(dir_name)
}

for (gene_info in Trichoblast_gene_list) {
  
  # create a feature plot for the gene
  feature_plot <- FeaturePlot(processed_LRP, features = gene_info[1], order = T, label = T, pt.size = 3) +
    scale_color_gradientn(colors = c('lightgray','yellow','red','darkred')) +
    ggtitle(paste0(gene_info[1], " - ", gene_info[2]))
  
  # Add the gene of interest to the metadata
  processed_LRP[["gene_of_interest"]] <- ifelse(processed_LRP[["RNA"]]@counts[gene_info[1],] > 0, gene_info[1], "Other")
  
  # create a dim plot for the gene
  dim_plot <- DimPlot(processed_LRP, group.by = "gene_of_interest") +
    ggtitle(paste0(gene_info[1], " - ", gene_info[2]))
  
  # Combine the plots into a single grob
  combined_plot <- gridExtra::arrangeGrob(feature_plot, dim_plot, ncol=2)
  
  # Specify the name of the PNG file
  png_file_name <- file.path(dir_name, paste0("Trichoblast_genes_", gene_info[1], "_", gene_info[2], ".png"))
  
  # Save the combined plot as a PNG file, specifying a larger width and height
  ggsave(filename = png_file_name, plot = combined_plot, width = 15, height = 5)
  
  # Remove the gene of interest from the metadata to prevent confusion in the next iteration
  processed_LRP[["gene_of_interest"]] <- NULL
}


# Define gene list
Epidermis_gene_list <- list(
  c("AT5G14750", "WER"),
  c("AT2G46410", "CPC")
)

# Define directory name
dir_name <- "Epidermis"

# Create directory if it doesn't exist
if(!dir.exists(dir_name)){
  dir.create(dir_name)
}

for (gene_info in Epidermis_gene_list) {
  
  # create a feature plot for the gene
  feature_plot <- FeaturePlot(processed_LRP, features = gene_info[1], order = T, label = T, pt.size = 3) +
    scale_color_gradientn(colors = c('lightgray','yellow','red','darkred')) +
    ggtitle(paste0(gene_info[1], " - ", gene_info[2]))
  
  # Add the gene of interest to the metadata
  processed_LRP[["gene_of_interest"]] <- ifelse(processed_LRP[["RNA"]]@counts[gene_info[1],] > 0, gene_info[1], "Other")
  
  # create a dim plot for the gene
  dim_plot <- DimPlot(processed_LRP, group.by = "gene_of_interest") +
    ggtitle(paste0(gene_info[1], " - ", gene_info[2]))
  
  # Combine the plots into a single grob
  combined_plot <- gridExtra::arrangeGrob(feature_plot, dim_plot, ncol=2)
  
  # Specify the name of the PNG file
  png_file_name <- file.path(dir_name, paste0("Epidermis_genes_", gene_info[1], "_", gene_info[2], ".png"))
  
  # Save the combined plot as a PNG file, specifying a larger width and height
  ggsave(filename = png_file_name, plot = combined_plot, width = 15, height = 5)
  
  # Remove the gene of interest from the metadata to prevent confusion in the next iteration
  processed_LRP[["gene_of_interest"]] <- NULL
}

# Define gene list
Cortex_gene_list <- list(
  c("AT1G05570", "CALS1"),
  c("AT1G15460", "BOR4"),
  c("AT1G62510", "AT1G62510"),
  c("AT3G12700", "AT3G12700"),
  c("AT3G14850", "TBL41"),
  c("AT5G03570", "IREG2"),
  c("AT5G44130", "FLA13"),
  c("AT5G53370", "PME61"),
  c("AT5G55250", "IAMT1"),
  c("AT1G64380", "ERF061"),
  c("AT5G48000", "CYP708A2"),
  c("AT3G21670", "NPF6.4"),
  c("AT1G30380", "PSAK"),
  c("AT1G52410", "TSA1"),
  c("AT3G26300", "CYP71B34"),
  c("AT4G00700", "AT4G00700"),
  c("AT3G56080", "AT3G56080"),
  c("AT5G27350", "SFP1"),
  c("AT2G32540", "CSLB4"),
  c("AT1G62500", "AT1G62500"),
  c("AT1G09750", "AED3"),
  c("AT1G62500", "CO2"),
  c("AT5G07990", "CYP75B1")
)

# Define directory name
dir_name <- "Cortex"

# Create directory if it doesn't exist
if(!dir.exists(dir_name)){
  dir.create(dir_name)
}

for (gene_info in Cortex_gene_list) {
  
  # create a feature plot for the gene
  feature_plot <- FeaturePlot(processed_LRP, features = gene_info[1], order = T, label = T, pt.size = 3) +
    scale_color_gradientn(colors = c('lightgray','yellow','red','darkred')) +
    ggtitle(paste0(gene_info[1], " - ", gene_info[2]))
  
  # Add the gene of interest to the metadata
  processed_LRP[["gene_of_interest"]] <- ifelse(processed_LRP[["RNA"]]@counts[gene_info[1],] > 0, gene_info[1], "Other")
  
  # create a dim plot for the gene
  dim_plot <- DimPlot(processed_LRP, group.by = "gene_of_interest") +
    ggtitle(paste0(gene_info[1], " - ", gene_info[2]))
  
  # Combine the plots into a single grob
  combined_plot <- gridExtra::arrangeGrob(feature_plot, dim_plot, ncol=2)
  
  # Specify the name of the PNG file
  png_file_name <- file.path(dir_name, paste0("Cortex_genes_", gene_info[1], "_", gene_info[2], ".png"))
  
  # Save the combined plot as a PNG file, specifying a larger width and height
  ggsave(filename = png_file_name, plot = combined_plot, width = 15, height = 5)
  
  # Remove the gene of interest from the metadata to prevent confusion in the next iteration
  processed_LRP[["gene_of_interest"]] <- NULL
}

# Define Endodermis gene list
Endodermis_gene_list <- list(
  c("AT1G30750", "AT1G30750"),
  c("AT1G60690", "AT1G60690"),
  c("AT1G61590", "AT1G61590"),
  c("AT1G71740", "AT1G71740"),
  c("AT2G21100", "DIR23"),
  c("AT2G27370", "CASP3"),
  c("AT2G30210", "LAC3"),
  c("AT2G39430", "DIR9"),
  c("AT2G40160", "TBL30"),
  c("AT2G48130", "AT2G48130"),
  c("AT3G11550", "CASP2"),
  c("AT3G22620", "AT3G22620"),
  c("AT4G02090", "AT4G02090"),
  c("AT4G17215", "AT4G17215"),
  c("AT5G06200", "CASP4"),
  c("AT5G42180", "PER64"),
  c("AT5G66390", "PER72"),
  c("AT1G01120", "KCS1"),
  c("AT2G14900", "GASA7"),
  c("AT1G07730", "DIR25"),
  c("AT2G48130", "AT2G48130"),
  c("AT3G45260", "IDD9"),
  c("AT1G61590", "PBL15"),
  c("AT5G57620", "MYB36"),
  c("AT2G31730", "BHLH154"),
  c("AT5G06200", "CASP4"),
  c("AT3G54220", "SCR"),
  c("AT4G28100", "EN7"),
  c("AT5G45210", "AT5G45210")
)
# Define directory name
dir_name <- "Endodermis"

# Create directory if it doesn't exist
if(!dir.exists(dir_name)){
  dir.create(dir_name)
}

for (gene_info in Endodermis_gene_list) {
  
  # create a feature plot for the gene
  feature_plot <- FeaturePlot(processed_LRP, features = gene_info[1], order = T, label = T, pt.size = 3) +
    scale_color_gradientn(colors = c('lightgray','yellow','red','darkred')) +
    ggtitle(paste0(gene_info[1], " - ", gene_info[2]))
  
  # Add the gene of interest to the metadata
  processed_LRP[["gene_of_interest"]] <- ifelse(processed_LRP[["RNA"]]@counts[gene_info[1],] > 0, gene_info[1], "Other")
  
  # create a dim plot for the gene
  dim_plot <- DimPlot(processed_LRP, group.by = "gene_of_interest") +
    ggtitle(paste0(gene_info[1], " - ", gene_info[2]))
  
  # Combine the plots into a single grob
  combined_plot <- gridExtra::arrangeGrob(feature_plot, dim_plot, ncol=2)
  
  # Specify the name of the PNG file
  png_file_name <- file.path(dir_name, paste0("Endodermis_genes_", gene_info[1], "_", gene_info[2], ".png"))
  
  # Save the combined plot as a PNG file, specifying a larger width and height
  ggsave(filename = png_file_name, plot = combined_plot, width = 15, height = 5)
  
  # Remove the gene of interest from the metadata to prevent confusion in the next iteration
  processed_LRP[["gene_of_interest"]] <- NULL
}

Pericycle_gene_list <- list(
  c("AT1G74430", "MYB95"),
  c("AT3G24310", "MYB305"),
  c("AT4G28250", "EXPB3"),
  c("AT5G25920", "AT5G25920"),
  c("AT4G29770", "AT4G29770"),
  c("AT1G58025", "AT1G58025"),
  c("AT4G34650", "SQS2"),
  c("AT5G06080", "LBD33"),
  c("AT1G52400", "BGLU18"),
  c("AT2G19900", "NADP-ME1"),
  c("AT5G24240", "PI4KG3"),
  c("AT3G06020", "FAF4"),
  c("AT3G29635", "AT3G29635"),
  c("AT4G17760", "AT4G17760"),
  c("AT4G28250", "EXPB3"),
  c("AT3G47320", "AT3G47320"),
  c("AT1G32450", "NPF7.3"),
  c("AT5G01740", "AT5G01740"),
  c("AT5G65420", "CYCD4;1"),
  c("AT3G23430", "PHO1"),
  c("AT5G26930", "GATA23"),
  c("AT1G11330", "RDA2")
)

# Define directory name
dir_name <- "Pericycle"

# Create directory if it doesn't exist
if(!dir.exists(dir_name)){
  dir.create(dir_name)
}

# loop through gene_ids to create a feature plot for each gene
for (gene_info in Pericycle_gene_list) {
  
  # create a feature plot for the gene
  feature_plot <- FeaturePlot(processed_LRP, features = gene_info[1], order = T, label = T, pt.size = 3) +
    scale_color_gradientn(colors = c('lightgray','yellow','red','darkred')) +
    ggtitle(paste0(gene_info[1], " - ", gene_info[2]))
  
  # Add the gene of interest to the metadata
  processed_LRP[["gene_of_interest"]] <- ifelse(processed_LRP[["RNA"]]@counts[gene_info[1],] > 0, gene_info[1], "Other")
  
  # create a dim plot for the gene
  dim_plot <- DimPlot(processed_LRP, group.by = "gene_of_interest") +
    ggtitle(paste0(gene_info[1], " - ", gene_info[2]))
  
  # Combine the plots into a single grob
  combined_plot <- gridExtra::arrangeGrob(feature_plot, dim_plot, ncol=2)
  
  # Specify the name of the PNG file
  png_file_name <- file.path(dir_name, paste0("Pericycle_genes_", gene_info[1], "_", gene_info[2], ".png"))
  
  # Save the combined plot as a PNG file, specifying a larger width and height
  ggsave(filename = png_file_name, plot = combined_plot, width = 15, height = 5)
  
  # Remove the gene of interest from the metadata to prevent confusion in the next iteration
  processed_LRP[["gene_of_interest"]] <- NULL
}

# Define Stele gene list
Stele_gene_list <- list(
  c("AT5G48070", "XTH20"),
  c("AT1G22710", "SUC2"),
  c("AT4G36710", "HAM4"),
  c("AT4G23800", "3xHGM-box2"),
  c("AT5G37660", "PDLP7"),
  c("AT4G21310", "AT4G21310"),
  c("AT4G36710", "SCL15"),
  c("AT4G37650", "SHR")
)

# Define directory name
dir_name <- "Stele"

# Create directory if it doesn't exist
if(!dir.exists(dir_name)){
  dir.create(dir_name)
}

# Loop through gene_ids to create a feature plot for each gene
for (gene_info in Stele_gene_list) {
  
  # create a feature plot for the gene
  feature_plot <- FeaturePlot(processed_LRP, features = gene_info[1], order = T, label = T, pt.size = 3) +
    scale_color_gradientn(colors = c('lightgray','yellow','red','darkred')) +
    ggtitle(paste0(gene_info[1], " - ", gene_info[2]))
  
  # Add the gene of interest to the metadata
  processed_LRP[["gene_of_interest"]] <- ifelse(processed_LRP[["RNA"]]@counts[gene_info[1],] > 0, gene_info[1], "Other")
  
  # create a dim plot for the gene
  dim_plot <- DimPlot(processed_LRP, group.by = "gene_of_interest") +
    ggtitle(paste0(gene_info[1], " - ", gene_info[2]))
  
  # Combine the plots into a single grob
  combined_plot <- gridExtra::arrangeGrob(feature_plot, dim_plot, ncol=2)
  
  # Specify the name of the PNG file
  png_file_name <- file.path(dir_name, paste0("Stele_genes_", gene_info[1], "_", gene_info[2], ".png"))
  
  # Save the combined plot as a PNG file, specifying a larger width and height
  ggsave(filename = png_file_name, plot = combined_plot, width = 15, height = 5)
  
  # Remove the gene of interest from the metadata to prevent confusion in the next iteration
  processed_LRP[["gene_of_interest"]] <- NULL
}

# Define gene list
Procambium_gene_list <- list(
  c("AT5G57130", "SMXL5"),
  c("AT4G32880", "ATHB-8"),
  c("AT3G09070", "OPS"),
  c("AT5G43810", "AGO10")
)

# Define directory name
dir_name <- "Procambium"

# Create directory if it doesn't exist
if(!dir.exists(dir_name)){
  dir.create(dir_name)
}

# Loop through gene list to create a feature plot for each gene
for (gene_info in Procambium_gene_list) {
  
  # create a feature plot for the gene
  feature_plot <- FeaturePlot(processed_LRP, features = gene_info[1], order = T, label = T, pt.size = 3) +
    scale_color_gradientn(colors = c('lightgray','yellow','red','darkred')) +
    ggtitle(paste0(gene_info[1], " - ", gene_info[2]))
  
  # Add the gene of interest to the metadata
  processed_LRP[["gene_of_interest"]] <- ifelse(processed_LRP[["RNA"]]@counts[gene_info[1],] > 0, gene_info[1], "Other")
  
  # create a dim plot for the gene
  dim_plot <- DimPlot(processed_LRP, group.by = "gene_of_interest") +
    ggtitle(paste0(gene_info[1], " - ", gene_info[2]))
  
  # Combine the plots into a single grob
  combined_plot <- gridExtra::arrangeGrob(feature_plot, dim_plot, ncol=2)
  
  # Specify the name of the PNG file
  png_file_name <- file.path(dir_name, paste0("Procambium_genes_", gene_info[1], "_", gene_info[2], ".png"))
  
  # Save the combined plot as a PNG file, specifying a larger width and height
  ggsave(filename = png_file_name, plot = combined_plot, width = 15, height = 5)
  
  # Remove the gene of interest from the metadata to prevent confusion in the next iteration
  processed_LRP[["gene_of_interest"]] <- NULL
}

# Define Metaxylem gene list
Metaxylem_gene_list <- list(
  c("AT1G80100", "AHP6"),
  c("AT4G37650", "SHR"),
  c("AT1G29950", "BHLH144"),
  c("AT4G02750", "PCMP-H24"),
  c("AT4G12620", "ORC1B"),
  c("AT1G74150", "AT1G74150"),
  c("AT5G08260", "SCPL35"),
  c("AT5G63140", "PAP29"),
  c("AT2G25060", "ENODL14"),
  c("AT3G14730", "PCMP-E31"),
  c("AT3G22790", "NET1A"),
  c("AT1G29940", "NRPA2"),
  c("AT5G06210", "AT5G06210"),
  c("AT3G61100", "AT3G61100"),
  c("AT2G13540", "ABH1"),
  c("AT3G53190", "AT3G53190"),
  c("AT5G37010", "AT5G37010"),
  c("AT5G25475", "AT5G25475"),
  c("AT3G54090", "FLN1"),
  c("AT3G20430", "AT3G20430"),
  c("AT3G25710", "TMO5")
)

# Define directory name
dir_name <- "Metaxylem"

# Create directory if it doesn't exist
if(!dir.exists(dir_name)){
  dir.create(dir_name)
}

# Loop through gene list to create a feature plot for each gene
for (gene_info in Metaxylem_gene_list) {
  
  # create a feature plot for the gene
  feature_plot <- FeaturePlot(processed_LRP, features = gene_info[1], order = T, label = T, pt.size = 3) +
    scale_color_gradientn(colors = c('lightgray','yellow','red','darkred')) +
    ggtitle(paste0(gene_info[1], " - ", gene_info[2]))
  
  # Add the gene of interest to the metadata
  processed_LRP[["gene_of_interest"]] <- ifelse(processed_LRP[["RNA"]]@counts[gene_info[1],] > 0, gene_info[1], "Other")
  
  # create a dim plot for the gene
  dim_plot <- DimPlot(processed_LRP, group.by = "gene_of_interest") +
    ggtitle(paste0(gene_info[1], " - ", gene_info[2]))
  
  # Combine the plots into a single grob
  combined_plot <- gridExtra::arrangeGrob(feature_plot, dim_plot, ncol=2)
  
  # Specify the name of the PNG file
  png_file_name <- file.path(dir_name, paste0("Metaxylem_genes_", gene_info[1], "_", gene_info[2], ".png"))
  
  # Save the combined plot as a PNG file, specifying a larger width and height
  ggsave(filename = png_file_name, plot = combined_plot, width = 15, height = 5)
  
  # Remove the gene of interest from the metadata to prevent confusion in the next iteration
  processed_LRP[["gene_of_interest"]] <- NULL
}

# Define gene list
Protoxylem_gene_list <- list(
  c("AT5G47635", "AT5G47635"),
  c("AT1G80100", "AHP6"),
  c("AT1G71930", "VND7"),
  c("AT1G04540", "AT1G04540"),
  c("AT1G04760", "VAMP726"),
  c("AT1G06020", "AT1G06020"),
  c("AT1G07120", "AT1G07120"),
  c("AT1G09610", "GXM1"),
  c("AT1G09890", "AT1G09890"),
  c("AT1G10800", "AT1G10800"),
  c("AT1G12450", "AT1G12450"),
  c("AT1G14190", "AT1G14190"),
  c("AT1G16520", "AT1G16520"),
  c("AT1G17950", "ATMYB52"),
  c("AT1G19230", "RBOHE"),
  c("AT1G20720", "AT1G20720"),
  c("AT1G21780", "AT1G21780"),
  c("AT1G26570", "UGD1"),
  c("AT1G27920", "MAP65-8"),
  c("AT1G29200", "AT1G29200"),
  c("AT1G48280", "AT1G48280"),
  c("AT1G51640", "ATEXO70G2"),
  c("AT1G66810", "AT1G66810"),
  c("AT5G62380", "NAC101"),
  c("AT1G68600", "ALMT5"),
  c("AT5G46115", "AT5G46115"),
  c("AT5G08480", "VQ31"),
  c("AT1G73410", "ATMYB54"),
  c("AT1G63910", "AtMYB103"),
  c("AT3G26125", "CYP86C2"),
  c("AT5G09480", "AT5G09480"),
  c("AT5G01730", "SCAR4"),
  c("AT3G52900", "AT3G52900"),
  c("AT5G11540", "GULLO3"),
  c("AT3G15050", "IQD10"),
  c("AT4G36540", "BEE2"),
  c("AT1G59850", "AT1G59850"),
  c("AT5G38610", "AT5G38610"),
  c("AT5G43070", "WPP1")
)

# Define directory name
dir_name <- "Protoxylem"

# Create directory if it doesn't exist
if(!dir.exists(dir_name)){
  dir.create(dir_name)
}

# loop through gene_ids to create a feature plot for each gene
for (gene_info in Protoxylem_gene_list ) {
  
  # create a feature plot for the gene
  feature_plot <- FeaturePlot(processed_LRP, features = gene_info[1], order = T, label = T, pt.size = 3) +
    scale_color_gradientn(colors = c('lightgray','yellow','red','darkred')) +
    ggtitle(paste0(gene_info[1], " - ", gene_info[2]))
  
  # Add the gene of interest to the metadata
  processed_LRP[["gene_of_interest"]] <- ifelse(processed_LRP[["RNA"]]@counts[gene_info[1],] > 0, gene_info[1], "Other")
  
  # create a dim plot for the gene
  dim_plot <- DimPlot(processed_LRP, group.by = "gene_of_interest") +
    ggtitle(paste0(gene_info[1], " - ", gene_info[2]))
  
  # Combine the plots into a single grob
  combined_plot <- gridExtra::arrangeGrob(feature_plot, dim_plot, ncol=2)
  
  # Specify the name of the PNG file
  png_file_name <- file.path(dir_name, paste0("Protoxylem_genes_", gene_info[1], "_", gene_info[2], ".png"))
  
  # Save the combined plot as a PNG file, specifying a larger width and height
  ggsave(filename = png_file_name, plot = combined_plot, width = 15, height = 5)
  
  # Remove the gene of interest from the metadata to prevent confusion in the next iteration
  processed_LRP[["gene_of_interest"]] <- NULL
}

# Define gene list for Protophloem
Protophloem_gene_list <- list(
  c("AT1G14730", "CYB561C"),
  c("AT5G40260", "SWEET8"),
  c("AT5G57130", "AT5G57130"),
  c("AT5G62940", "DOF5.6"),
  c("AT1G15080", "LPP2"),
  c("AT3G01680", "SEOB"),
  c("AT4G13600", "AT4G13600"),
  c("AT5G67460", "AT5G67460"),
  c("AT5G22870", "AT5G22870"),
  c("AT1G54330", "ANAC020"),
  c("AT1G60030", "NAT7"),
  c("AT3G17730", "anac057"),
  c("AT4G17910", "AT4G17910"),
  c("AT4G36620", "GATA19"),
  c("AT4G03270", "CYCD6-1"),
  c("AT1G06490", "CALS7"),
  c("AT1G11915", "AT1G11915"),
  c("AT1G13050", "AT1G13050"),
  c("AT1G31335", "AT1G31335"),
  c("AT1G33480", "ATL58"),
  c("AT5G01370", "ACI1")
)

# Define directory name
dir_name <- "Protophloem"

# Create directory if it doesn't exist
if(!dir.exists(dir_name)){
  dir.create(dir_name)
}

# loop through gene_ids to create a feature plot for each gene
for (gene_info in Protophloem_gene_list ) {
  
  # create a feature plot for the gene
  feature_plot <- FeaturePlot(processed_LRP, features = gene_info[1], order = T, label = T, pt.size = 3) +
    scale_color_gradientn(colors = c('lightgray','yellow','red','darkred')) +
    ggtitle(paste0(gene_info[1], " - ", gene_info[2]))
  
  # Add the gene of interest to the metadata
  processed_LRP[["gene_of_interest"]] <- ifelse(processed_LRP[["RNA"]]@counts[gene_info[1],] > 0, gene_info[1], "Other")
  
  # create a dim plot for the gene
  dim_plot <- DimPlot(processed_LRP, group.by = "gene_of_interest") +
    ggtitle(paste0(gene_info[1], " - ", gene_info[2]))
  
  # Combine the plots into a single grob
  combined_plot <- gridExtra::arrangeGrob(feature_plot, dim_plot, ncol=2)
  
  # Specify the name of the PNG file
  png_file_name <- file.path(dir_name, paste0("Protophloem_genes_", gene_info[1], "_", gene_info[2], ".png"))
  
  # Save the combined plot as a PNG file, specifying a larger width and height
  ggsave(filename = png_file_name, plot = combined_plot, width = 15, height = 5)
  
  # Remove the gene of interest from the metadata to prevent confusion in the next iteration
  processed_LRP[["gene_of_interest"]] <- NULL
}

# Define the Phloem gene list
Phloem_gene_list <- list(
  c("AT1G11570", "NTL"),
  c("AT1G69580", "AT1G69580"),
  c("AT1G69970", "CLE26"),
  c("AT2G41200", "AT2G41200"),
  c("AT3G47180", "AT3G47180"),
  c("AT4G04620", "ATG8B"),
  c("AT4G16400", "AT4G16400"),
  c("AT5G01480", "AT5G01480"),
  c("AT5G04680", "AT5G04680"),
  c("AT5G16610", "AT5G16610"),
  c("AT5G53030", "AT5G53030"),
  c("AT1G02950", "GSTF4"),
  c("AT1G05610", "APS2"),
  c("AT1G05770", "JAL2"),
  c("AT1G10360", "GSTU18"),
  c("AT1G13920", "AT1G13920"),
  c("AT5G20030", "AT5G20030"),
  c("AT3G58350", "RTM3"),
  c("AT1G79430", "APL"),
  c("AT5G62940", "DOF5.6"),
  c("AT3G43270", "PME32"),
  c("AT1G26790", "CDF 6"),
  c("AT1G29520", "AT1G29520"),
  c("AT2G37590", "DOF2.4"),
  c("AT1G07640", "OBP2")
)

# Define directory name
dir_name <- "Phloem"

# Create directory if it doesn't exist
if(!dir.exists(dir_name)){
  dir.create(dir_name)
}

# loop through gene_ids to create a feature plot for each gene
for (gene_info in Phloem_gene_list ) {
  
  # create a feature plot for the gene
  feature_plot <- FeaturePlot(processed_LRP, features = gene_info[1], order = T, label = T, pt.size = 3) +
    scale_color_gradientn(colors = c('lightgray','yellow','red','darkred')) +
    ggtitle(paste0(gene_info[1], " - ", gene_info[2]))
  
  # Add the gene of interest to the metadata
  processed_LRP[["gene_of_interest"]] <- ifelse(processed_LRP[["RNA"]]@counts[gene_info[1],] > 0, gene_info[1], "Other")
  
  # create a dim plot for the gene
  dim_plot <- DimPlot(processed_LRP, group.by = "gene_of_interest") +
    ggtitle(paste0(gene_info[1], " - ", gene_info[2]))
  
  # Combine the plots into a single grob
  combined_plot <- gridExtra::arrangeGrob(feature_plot, dim_plot, ncol=2)
  
  # Specify the name of the PNG file
  png_file_name <- file.path(dir_name, paste0("Phloem_genes_", gene_info[1], "_", gene_info[2], ".png"))
  
  # Save the combined plot as a PNG file, specifying a larger width and height
  ggsave(filename = png_file_name, plot = combined_plot, width = 15, height = 5)
  
  # Remove the gene of interest from the metadata to prevent confusion in the next iteration
  processed_LRP[["gene_of_interest"]] <- NULL
}

# Define gene list
PSE_gene_list <- list(
  c("AT5G17260", "NAC086"),
  c("AT3G03200", "NAC45"),
  c("AT3G01680", "ATSEOR1"),
  c("AT2G37590", "DOF2.4"),
  c("AT3G22400", "LOX5")
)

# Define directory name
dir_name <- "PSE"

# Create directory if it doesn't exist
if(!dir.exists(dir_name)){
  dir.create(dir_name)
}

# loop through gene_ids to create a feature plot for each gene
for (gene_info in PSE_gene_list ) {
  
  # create a feature plot for the gene
  feature_plot <- FeaturePlot(processed_LRP, features = gene_info[1], order = T, label = T, pt.size = 3) +
    scale_color_gradientn(colors = c('lightgray','yellow','red','darkred')) +
    ggtitle(paste0(gene_info[1], " - ", gene_info[2]))
  
  # Add the gene of interest to the metadata
  processed_LRP[["gene_of_interest"]] <- ifelse(processed_LRP[["RNA"]]@counts[gene_info[1],] > 0, gene_info[1], "Other")
  
  # create a dim plot for the gene
  dim_plot <- DimPlot(processed_LRP, group.by = "gene_of_interest") +
    ggtitle(paste0(gene_info[1], " - ", gene_info[2]))
  
  # Combine the plots into a single grob
  combined_plot <- gridExtra::arrangeGrob(feature_plot, dim_plot, ncol=2)
  
  # Specify the name of the PNG file
  png_file_name <- file.path(dir_name, paste0("PSE_genes_", gene_info[1], "_", gene_info[2], ".png"))
  
  # Save the combined plot as a PNG file, specifying a larger width and height
  ggsave(filename = png_file_name, plot = combined_plot, width = 15, height = 5)
  
  # Remove the gene of interest from the metadata to prevent confusion in the next iteration
  processed_LRP[["gene_of_interest"]] <- NULL
}

# Define the PCC gene list
PCC_gene_list <- list(
  c("AT5G17260", "NAC086"),
  c("AT3G03200", "NAC45"),
  c("AT3G01680", "ATSEOR1"),
  c("AT2G37590", "DOF2.4"),
  c("AT3G22400", "LOX5")
)

# Define directory name
dir_name <- "PCC"

# Create directory if it doesn't exist
if(!dir.exists(dir_name)){
  dir.create(dir_name)
}

# loop through gene_ids to create a feature plot for each gene
for (gene_info in PCC_gene_list ) {
  
  # create a feature plot for the gene
  feature_plot <- FeaturePlot(processed_LRP, features = gene_info[1], order = T, label = T, pt.size = 3) +
    scale_color_gradientn(colors = c('lightgray','yellow','red','darkred')) +
    ggtitle(paste0(gene_info[1], " - ", gene_info[2]))
  
  # Add the gene of interest to the metadata
  processed_LRP[["gene_of_interest"]] <- ifelse(processed_LRP[["RNA"]]@counts[gene_info[1],] > 0, gene_info[1], "Other")
  
  # create a dim plot for the gene
  dim_plot <- DimPlot(processed_LRP, group.by = "gene_of_interest") +
    ggtitle(paste0(gene_info[1], " - ", gene_info[2]))
  
  # Combine the plots into a single grob
  combined_plot <- gridExtra::arrangeGrob(feature_plot, dim_plot, ncol=2)
  
  # Specify the name of the PNG file
  png_file_name <- file.path(dir_name, paste0("PCC_genes_", gene_info[1], "_", gene_info[2], ".png"))
  
  # Save the combined plot as a PNG file, specifying a larger width and height
  ggsave(filename = png_file_name, plot = combined_plot, width = 15, height = 5)
  
  # Remove the gene of interest from the metadata to prevent confusion in the next iteration
  processed_LRP[["gene_of_interest"]] <- NULL
}

# Define gene list
QC_gene_list <- list(
  c("AT1G68640", "PAN"),
  c("AT3G26120", "TEL1"),
  c("AT5G23780", "DUF9"),
  c("AT5G17430", "BBM"),
  c("AT2G03830", "RGF8"),
  c("AT2G23050", "NPY4"),
  c("AT1G30530", "UGT78D1"),
  c("AT5G25900", "KO"),
  c("AT5G45780", "AT5G45780"),
  c("AT1G11420", "ATDUF2"),
  c("AT3G13510", "AT3G13510"),
  c("AT3G30350", "RGF4"),
  c("AT1G16070", "TULP8"),
  c("AT1G75580", "AT1G75580"),
  c("AT3G11260", "WOX5"),
  c("AT5G62165", "AGL42"),
  c("AT1G55200", "AT1G55200")
)

# Define directory name
dir_name <- "QC"

# Create directory if it doesn't exist
if(!dir.exists(dir_name)){
  dir.create(dir_name)
}

# loop through gene_ids to create a feature plot for each gene
for (gene_info in QC_gene_list ) {
  
  # create a feature plot for the gene
  feature_plot <- FeaturePlot(processed_LRP, features = gene_info[1], order = T, label = T, pt.size = 3) +
    scale_color_gradientn(colors = c('lightgray','yellow','red','darkred')) +
    ggtitle(paste0(gene_info[1], " - ", gene_info[2]))
  
  # Add the gene of interest to the metadata
  processed_LRP[["gene_of_interest"]] <- ifelse(processed_LRP[["RNA"]]@counts[gene_info[1],] > 0, gene_info[1], "Other")
  
  # create a dim plot for the gene
  dim_plot <- DimPlot(processed_LRP, group.by = "gene_of_interest") +
    ggtitle(paste0(gene_info[1], " - ", gene_info[2]))
  
  # Combine the plots into a single grob
  combined_plot <- gridExtra::arrangeGrob(feature_plot, dim_plot, ncol=2)
  
  # Specify the name of the PNG file
  png_file_name <- file.path(dir_name, paste0("QC_genes_", gene_info[1], "_", gene_info[2], ".png"))
  
  # Save the combined plot as a PNG file, specifying a larger width and height
  ggsave(filename = png_file_name, plot = combined_plot, width = 15, height = 5)
  
  # Remove the gene of interest from the metadata to prevent confusion in the next iteration
  processed_LRP[["gene_of_interest"]] <- NULL
}

Pericycle_xylem_pole_gene_list <- list(
  c("AT2G18800", "XTH21"),
  c("AT5G52470", "FIB1"),
  c("AT4G18250", "AT4G18250"),
  c("AT5G01740", "AT5G01740"),
  c("AT5G26930", "GATA23")
)

# Define directory name
dir_name <- "Pericycle_xylem_pole_genes"

# Create directory if it doesn't exist
if(!dir.exists(dir_name)){
  dir.create(dir_name)
}

# loop through gene_ids to create a feature plot for each gene
for (gene_info in Pericycle_xylem_pole_gene_list ) {
  
  # create a feature plot for the gene
  plot <- FeaturePlot(processed_LRP, features = gene_info[1], order = T, label = T, pt.size = 3) +
    scale_color_gradientn(colors = c('lightgray','yellow','red','darkred')) +
    ggtitle(paste0(gene_info[1], " - ", gene_info[2]))
  
  # Specify the name of the PDF file
  pdf_file_name <- file.path(dir_name, paste0("Pericycle_xylem_pole_gene_list_", gene_info[1], "_", gene_info[2], ".pdf"))
  
  # Start PDF file
  pdf(file = pdf_file_name)
  
  # Print the plot to the file
  print(plot)
  
  # Close the PDF file
  dev.off()
}

Pericycle_phloem_pole_gene_list <- list(
  c("AT5G01750", "AT5G01750"),
  c("AT2G22850", "AT2G22850"),
  c("AT1G18570", "MYB51"),
  c("AT1G35910", "TPPD"),
  c("AT4G23410", "TET5")
)

# Define directory name
dir_name <- "Pericycle_phloem_pole"

# Create directory if it doesn't exist
if(!dir.exists(dir_name)){
  dir.create(dir_name)
}

# loop through gene_ids to create a feature plot for each gene
for (gene_info in Pericycle_phloem_pole_gene_list ) {
  
  # create a feature plot for the gene
  feature_plot <- FeaturePlot(processed_LRP, features = gene_info[1], order = T, label = T, pt.size = 3) +
    scale_color_gradientn(colors = c('lightgray','yellow','red','darkred')) +
    ggtitle(paste0(gene_info[1], " - ", gene_info[2]))
  
  # Add the gene of interest to the metadata
  processed_LRP[["gene_of_interest"]] <- ifelse(processed_LRP[["RNA"]]@counts[gene_info[1],] > 0, gene_info[1], "Other")
  
  # create a dim plot for the gene
  dim_plot <- DimPlot(processed_LRP, group.by = "gene_of_interest") +
    ggtitle(paste0(gene_info[1], " - ", gene_info[2]))
  
  # Combine the plots into a single grob
  combined_plot <- gridExtra::arrangeGrob(feature_plot, dim_plot, ncol=2)
  
  # Specify the name of the PNG file
  png_file_name <- file.path(dir_name, paste0("Pericycle_phloem_pole_genes_", gene_info[1], "_", gene_info[2], ".png"))
  
  # Save the combined plot as a PNG file, specifying a larger width and height
  ggsave(filename = png_file_name, plot = combined_plot, width = 15, height = 5)
  
  # Remove the gene of interest from the metadata to prevent confusion in the next iteration
  processed_LRP[["gene_of_interest"]] <- NULL
}

# Define gene list
LRC_gene_list <- list(
  c("AT1G33280", "BRN1"),
  c("AT4G10350", "BRN2"),
  c("AT1G79580", "SMB"),
  c("AT2G32610", "CSLB01"),
  c("AT4G37160", "sks15"),
  c("AT3G55930", "AT3G55930"),
  c("AT5G36880", "ACS"),
  c("AT1G54010", "GLL23"),
  c("AT1G14220", "AT1G14220")
)

# Define directory name
dir_name <- "LRC_gene_list"

# Create directory if it doesn't exist
if(!dir.exists(dir_name)){
  dir.create(dir_name)
}

# loop through gene_ids to create a feature plot for each gene
for (gene_info in LRC_gene_list ) {
  
  # create a feature plot for the gene
  feature_plot <- FeaturePlot(processed_LRP, features = gene_info[1], order = T, label = T, pt.size = 3) +
    scale_color_gradientn(colors = c('lightgray','yellow','red','darkred')) +
    ggtitle(paste0(gene_info[1], " - ", gene_info[2]))
  
  # Add the gene of interest to the metadata
  processed_LRP[["gene_of_interest"]] <- ifelse(processed_LRP[["RNA"]]@counts[gene_info[1],] > 0, gene_info[1], "Other")
  
  # create a dim plot for the gene
  dim_plot <- DimPlot(processed_LRP, group.by = "gene_of_interest") +
    ggtitle(paste0(gene_info[1], " - ", gene_info[2]))
  
  # Combine the plots into a single grob
  combined_plot <- gridExtra::arrangeGrob(feature_plot, dim_plot, ncol=2)
  
  # Specify the name of the PNG file
  png_file_name <- file.path(dir_name, paste0("LRC_genes_", gene_info[1], "_", gene_info[2], ".png"))
  
  # Save the combined plot as a PNG file, specifying a larger width and height
  ggsave(filename = png_file_name, plot = combined_plot, width = 15, height = 5)
  
  # Remove the gene of interest from the metadata to prevent confusion in the next iteration
  processed_LRP[["gene_of_interest"]] <- NULL
}

# Define gene list
Columella_gene_list <- list(
  c("AT2G33230", "YUC7"),
  c("AT3G10320", ""),
  c("AT4G18350", "NCED2"),
  c("AT5G22550", ""),
  c("AT1G26680", "REM17"),
  c("AT2G18230", "PPA2"),
  c("AT3G18230", ""),
  c("AT3G28540", ""),
  c("AT3G52890", "KIPK"),
  c("AT3G55180", ""),
  c("AT5G07550", "GRP19"),
  c("AT3G29810", "COBL3"),
  c("AT5G48300", "APS1"),
  c("AT5G64860", "DPE1"),
  c("AT3G44240", "CAF1-8"),
  c("AT3G62320", ""),
  c("AT5G48130", ""),
  c("AT1G04600", "XI-A"),
  c("AT5G58580", "ATL63"),
  c("AT4G09960", "AGL11"),
  c("AT1G69830", "AMY3"),
  c("AT3G28540", ""),
  c("AT1G26870", "FEZ"),
  c("AT3G62320", ""),
  c("AT3G10600", "CAT7"),
  c("AT1G17400", "LZY2"),
  c("AT1G19115", "LZY3")
)

# Define directory name
dir_name <- "Columella_gene_list"

# Create directory if it doesn't exist
if(!dir.exists(dir_name)){
  dir.create(dir_name)
}

# loop through gene_ids to create a feature plot for each gene
for (gene_info in Columella_gene_list ) {
  
  # create a feature plot for the gene
  feature_plot <- FeaturePlot(processed_LRP, features = gene_info[1], order = T, label = T, pt.size = 3) +
    scale_color_gradientn(colors = c('lightgray','yellow','red','darkred')) +
    ggtitle(paste0(gene_info[1], " - ", gene_info[2]))
  
  # Add the gene of interest to the metadata
  processed_LRP[["gene_of_interest"]] <- ifelse(processed_LRP[["RNA"]]@counts[gene_info[1],] > 0, gene_info[1], "Other")
  
  # create a dim plot for the gene
  dim_plot <- DimPlot(processed_LRP, group.by = "gene_of_interest") +
    ggtitle(paste0(gene_info[1], " - ", gene_info[2]))
  
  # Combine the plots into a single grob
  combined_plot <- gridExtra::arrangeGrob(feature_plot, dim_plot, ncol=2)
  
  # Specify the name of the PNG file
  png_file_name <- file.path(dir_name, paste0("Columella_genes_", gene_info[1], "_", gene_info[2], ".png"))
  
  # Save the combined plot as a PNG file, specifying a larger width and height
  ggsave(filename = png_file_name, plot = combined_plot, width = 15, height = 5)
  
  # Remove the gene of interest from the metadata to prevent confusion in the next iteration
  processed_LRP[["gene_of_interest"]] <- NULL
}

# Define gene list
stele_mapping_gene_list <- list(
  c("AT3G48100", "ARR5"),
  c("AT3G23030", "IAA2"),
  c("AT1G80100", "AHP6"),
  c("AT4G32880", "ATHB8"),
  c("AT1G73590", "PIN1"),
  c("AT1G70940", "PIN3"),
  c("AT1G23080", "PIN7"),
  c("AT1G71930", "VND7"),
  c("AT2G38120", "Aux1"),
  c("AT2G34710", "PHB")
)

# Define directory name
dir_name <- "stele_mapping_gene"

# Create directory if it doesn't exist
if(!dir.exists(dir_name)){
  dir.create(dir_name)
}

# loop through gene_ids to create a feature plot for each gene
for (gene_info in stele_mapping_gene_list ) {
  
  # create a feature plot for the gene
  feature_plot <- FeaturePlot(processed_LRP, features = gene_info[1], order = T, label = T, pt.size = 3) +
    scale_color_gradientn(colors = c('lightgray','yellow','red','darkred')) +
    ggtitle(paste0(gene_info[1], " - ", gene_info[2]))
  
  # Add the gene of interest to the metadata
  processed_LRP[["gene_of_interest"]] <- ifelse(processed_LRP[["RNA"]]@counts[gene_info[1],] > 0, gene_info[1], "Other")
  
  # create a dim plot for the gene
  dim_plot <- DimPlot(processed_LRP, group.by = "gene_of_interest") +
    ggtitle(paste0(gene_info[1], " - ", gene_info[2]))
  
  # Combine the plots into a single grob
  combined_plot <- gridExtra::arrangeGrob(feature_plot, dim_plot, ncol=2)
  
  # Specify the name of the PNG file
  png_file_name <- file.path(dir_name, paste0("Stele_mapping_genes_", gene_info[1], "_", gene_info[2], ".png"))
  
  # Save the combined plot as a PNG file, specifying a larger width and height
  ggsave(filename = png_file_name, plot = combined_plot, width = 15, height = 5)
  
  # Remove the gene of interest from the metadata to prevent confusion in the next iteration
  processed_LRP[["gene_of_interest"]] <- NULL
}

# Define gene list
Initials_gene_list <- list(
  c("AT3G60390", "HAT3"),
  c("AT3G56070", "CYP19-3"),
  c("AT5G10400", "H3.1"),
  c("AT2G27250", "CLV3")
)

# Define directory name
dir_name <- "Initials_genes"

# Create directory if it doesn't exist
if(!dir.exists(dir_name)){
  dir.create(dir_name)
}

# loop through gene_ids to create a feature plot for each gene
for (gene_info in Initials_gene_list ) {
  
  # create a feature plot for the gene
  feature_plot <- FeaturePlot(processed_LRP, features = gene_info[1], order = T, label = T, pt.size = 3) +
    scale_color_gradientn(colors = c('lightgray','yellow','red','darkred')) +
    ggtitle(paste0(gene_info[1], " - ", gene_info[2]))
  
  # Add the gene of interest to the metadata
  processed_LRP[["gene_of_interest"]] <- ifelse(processed_LRP[["RNA"]]@counts[gene_info[1],] > 0, gene_info[1], "Other")
  
  # create a dim plot for the gene
  dim_plot <- DimPlot(processed_LRP, group.by = "gene_of_interest") +
    ggtitle(paste0(gene_info[1], " - ", gene_info[2]))
  
  # Combine the plots into a single grob
  combined_plot <- gridExtra::arrangeGrob(feature_plot, dim_plot, ncol=2)
  
  # Specify the name of the PNG file
  png_file_name <- file.path(dir_name, paste0("Initials_genes_", gene_info[1], "_", gene_info[2], ".png"))
  
  # Save the combined plot as a PNG file, specifying a larger width and height
  ggsave(filename = png_file_name, plot = combined_plot, width = 15, height = 5)
  
  # Remove the gene of interest from the metadata to prevent confusion in the next iteration
  processed_LRP[["gene_of_interest"]] <- NULL
}

# Define gene list
VT_gene_list <- list(
  c("AT2G01830", "WOL")
)

# Define directory name
dir_name <- "VT"

# Create directory if it doesn't exist
if(!dir.exists(dir_name)){
  dir.create(dir_name)
}

# Loop through gene_ids to create a feature plot for each gene
for (gene_info in VT_gene_list ) {
  
  # create a feature plot for the gene
  feature_plot <- FeaturePlot(processed_LRP, features = gene_info[1], order = T, label = T, pt.size = 3) +
    scale_color_gradientn(colors = c('lightgray','yellow','red','darkred')) +
    ggtitle(paste0(gene_info[1], " - ", gene_info[2]))
  
  # Add the gene of interest to the metadata
  processed_LRP[["gene_of_interest"]] <- ifelse(processed_LRP[["RNA"]]@counts[gene_info[1],] > 0, gene_info[1], "Other")
  
  # create a dim plot for the gene
  dim_plot <- DimPlot(processed_LRP, group.by = "gene_of_interest") +
    ggtitle(paste0(gene_info[1], " - ", gene_info[2]))
  
  # Combine the plots into a single grob
  combined_plot <- gridExtra::arrangeGrob(feature_plot, dim_plot, ncol=2)
  
  # Specify the name of the PNG file
  png_file_name <- file.path(dir_name, paste0("VT_genes_", gene_info[1], "_", gene_info[2], ".png"))
  
  # Save the combined plot as a PNG file, specifying a larger width and height
  ggsave(filename = png_file_name, plot = combined_plot, width = 15, height = 5)
  
  # Remove the gene of interest from the metadata to prevent confusion in the next iteration
  processed_LRP[["gene_of_interest"]] <- NULL
}
# Define gene list
Cell_D_gene_list <- list(
  c("AT1G11190", "ENDO1"),
  c("AT1G26820", "RNS3"),
  c("AT3G45010", "SCPL48"),
  c("AT4G18425", "DMP4"),
  c("AT5G50260", "CEP1"),
  c("AT2G14095", "EXI1")
)

# Define directory name
dir_name <- "Cell_D"

# Create directory if it doesn't exist
if(!dir.exists(dir_name)){
  dir.create(dir_name)
}

# loop through gene_ids to create a feature plot for each gene
for (gene_info in Cell_D_gene_list ) {
  
  # create a feature plot for the gene
  feature_plot <- FeaturePlot(processed_LRP, features = gene_info[1], order = T, label = T, pt.size = 3) +
    scale_color_gradientn(colors = c('lightgray','yellow','red','darkred')) +
    ggtitle(paste0(gene_info[1], " - ", gene_info[2]))
  
  # Add the gene of interest to the metadata
  processed_LRP[["gene_of_interest"]] <- ifelse(processed_LRP[["RNA"]]@counts[gene_info[1],] > 0, gene_info[1], "Other")
  
  # create a dim plot for the gene
  dim_plot <- DimPlot(processed_LRP, group.by = "gene_of_interest") +
    ggtitle(paste0(gene_info[1], " - ", gene_info[2]))
  
  # Combine the plots into a single grob
  combined_plot <- gridExtra::arrangeGrob(feature_plot, dim_plot, ncol=2)
  
  # Specify the name of the PNG file
  png_file_name <- file.path(dir_name, paste0("Cell_D_genes_", gene_info[1], "_", gene_info[2], ".png"))
  
  # Save the combined plot as a PNG file, specifying a larger width and height
  ggsave(filename = png_file_name, plot = combined_plot, width = 15, height = 5)
  
  # Remove the gene of interest from the metadata to prevent confusion in the next iteration
  processed_LRP[["gene_of_interest"]] <- NULL
}

# Define gene list
Metaphloem_gene_list <- list(
  c("AT1G61760", "SEMA2")
)

# Define directory name
dir_name <- "Metaphloem"

# Create directory if it doesn't exist
if(!dir.exists(dir_name)){
  dir.create(dir_name)
}

# loop through gene_ids to create a feature plot for each gene
for (gene_info in Metaphloem_gene_list ) {
  
  # create a feature plot for the gene
  feature_plot <- FeaturePlot(processed_LRP, features = gene_info[1], order = T, label = T, pt.size = 3) +
    scale_color_gradientn(colors = c('lightgray','yellow','red','darkred')) +
    ggtitle(paste0(gene_info[1], " - ", gene_info[2]))
  
  # Add the gene of interest to the metadata
  processed_LRP[["gene_of_interest"]] <- ifelse(processed_LRP[["RNA"]]@counts[gene_info[1],] > 0, gene_info[1], "Other")
  
  # create a dim plot for the gene
  dim_plot <- DimPlot(processed_LRP, group.by = "gene_of_interest") +
    ggtitle(paste0(gene_info[1], " - ", gene_info[2]))
  
  # Combine the plots into a single grob
  combined_plot <- gridExtra::arrangeGrob(feature_plot, dim_plot, ncol=2)
  
  # Specify the name of the PNG file
  png_file_name <- file.path(dir_name, paste0("Metaphloem_genes_", gene_info[1], "_", gene_info[2], ".png"))
  
  # Save the combined plot as a PNG file, specifying a larger width and height
  ggsave(filename = png_file_name, plot = combined_plot, width = 15, height = 5)
  
  # Remove the gene of interest from the metadata to prevent confusion in the next iteration
  processed_LRP[["gene_of_interest"]] <- NULL
}

# Define gene list
SCN_gene_list <- list(
  c("AT5G51451", "RGF5"),
  c("AT5G60810", "RGF1"),
  c("AT2G03830", "RGF8"),
  c("AT5G17800", "BRAVO"),
  c("AT3G20880","WIP4")
)

# Define directory name
dir_name <- "SCN"

# Create directory if it doesn't exist
if(!dir.exists(dir_name)){
  dir.create(dir_name)
}

# loop through gene_ids to create a feature plot for each gene
for (gene_info in SCN_gene_list ) {
  
  # create a feature plot for the gene
  feature_plot <- FeaturePlot(processed_LRP, features = gene_info[1], order = T, label = T, pt.size = 3) +
    scale_color_gradientn(colors = c('lightgray','yellow','red','darkred')) +
    ggtitle(paste0(gene_info[1], " - ", gene_info[2]))
  
  # Add the gene of interest to the metadata
  processed_LRP[["gene_of_interest"]] <- ifelse(processed_LRP[["RNA"]]@counts[gene_info[1],] > 0, gene_info[1], "Other")
  
  # create a dim plot for the gene
  dim_plot <- DimPlot(processed_LRP, group.by = "gene_of_interest") +
    ggtitle(paste0(gene_info[1], " - ", gene_info[2]))
  
  # Combine the plots into a single grob
  combined_plot <- gridExtra::arrangeGrob(feature_plot, dim_plot, ncol=2)
  
  # Specify the name of the PNG file
  png_file_name <- file.path(dir_name, paste0("SCN_genes_", gene_info[1], "_", gene_info[2], ".png"))
  
  # Save the combined plot as a PNG file, specifying a larger width and height
  ggsave(filename = png_file_name, plot = combined_plot, width = 15, height = 5)
  
  # Remove the gene of interest from the metadata to prevent confusion in the next iteration
  processed_LRP[["gene_of_interest"]] <- NULL
}

#Define gene list
Dividing_cells_gene_list <- list(
  c("AT1G08560", "KN"),
  c("AT1G76540", "CDKB2-1"),
  c("AT1G44110", "CYCA1-1"),
  c("AT4G37490", "CYCB1-1"),
  c("AT2G17620", "CYCB2-1")
)

# Define directory name
dir_name <- "Dividing_cells"

# Create directory if it doesn't exist
if(!dir.exists(dir_name)){
  dir.create(dir_name)
}

# loop through gene_ids to create a feature plot for each gene
for (gene_info in Dividing_cells_gene_list ) {
  
  # create a feature plot for the gene
  feature_plot <- FeaturePlot(processed_LRP, features = gene_info[1], order = T, label = T, pt.size = 3) +
    scale_color_gradientn(colors = c('lightgray','yellow','red','darkred')) +
    ggtitle(paste0(gene_info[1], " - ", gene_info[2]))
  
  # Add the gene of interest to the metadata
  processed_LRP[["gene_of_interest"]] <- ifelse(processed_LRP[["RNA"]]@counts[gene_info[1],] > 0, gene_info[1], "Other")
  
  # create a dim plot for the gene
  dim_plot <- DimPlot(processed_LRP, group.by = "gene_of_interest") +
    ggtitle(paste0(gene_info[1], " - ", gene_info[2]))
  
  # Combine the plots into a single grob
  combined_plot <- gridExtra::arrangeGrob(feature_plot, dim_plot, ncol=2)
  
  # Specify the name of the PNG file
  png_file_name <- file.path(dir_name, paste0("Dividing_cells_genes_", gene_info[1], "_", gene_info[2], ".png"))
  
  # Save the combined plot as a PNG file, specifying a larger width and height
  ggsave(filename = png_file_name, plot = combined_plot, width = 15, height = 5)
  
  # Remove the gene of interest from the metadata to prevent confusion in the next iteration
  processed_LRP[["gene_of_interest"]] <- NULL
}

# Define gene list
Zone_markers_m_gene_list <- list(
  c("AT4G28190", "ULT1"),
  c("AT1G04020", "ROW1-BARD1"),
  c("AT4G36930", "SPT"),
  c("AT1G74500", "TMO7"),
  c("AT3G19380", "PUB25")
)

# Define directory name
dir_name <- "Zone_markers_m"

# Create directory if it doesn't exist
if(!dir.exists(dir_name)){
  dir.create(dir_name)
}

# Loop through gene_ids to create a feature plot for each gene
for (gene_info in Zone_markers_m_gene_list ) {
  
  # create a feature plot for the gene
  feature_plot <- FeaturePlot(processed_LRP, features = gene_info[1], order = T, label = T, pt.size = 3) +
    scale_color_gradientn(colors = c('lightgray','yellow','red','darkred')) +
    ggtitle(paste0(gene_info[1], " - ", gene_info[2]))
  
  # Add the gene of interest to the metadata
  processed_LRP[["gene_of_interest"]] <- ifelse(processed_LRP[["RNA"]]@counts[gene_info[1],] > 0, gene_info[1], "Other")
  
  # create a dim plot for the gene
  dim_plot <- DimPlot(processed_LRP, group.by = "gene_of_interest") +
    ggtitle(paste0(gene_info[1], " - ", gene_info[2]))
  
  # Combine the plots into a single grob
  combined_plot <- gridExtra::arrangeGrob(feature_plot, dim_plot, ncol=2)
  
  # Specify the name of the PNG file
  png_file_name <- file.path(dir_name, paste0("Zone_markers_M_genes_", gene_info[1], "_", gene_info[2], ".png"))
  
  # Save the combined plot as a PNG file, specifying a larger width and height
  ggsave(filename = png_file_name, plot = combined_plot, width = 15, height = 5)
  
  # Remove the gene of interest from the metadata to prevent confusion in the next iteration
  processed_LRP[["gene_of_interest"]] <- NULL
}

# Define gene list
Zone_markers_t_gene_list <- list(
  c("AT4G32040", "KNAT5"),
  c("AT5G04470", "SIM"),
  c("AT3G10525", "SMR1"),
  c("AT4G22910", "CCS52A1")
)

# Define directory name
dir_name <- "Zone_markers_t"

# Create directory if it doesn't exist
if(!dir.exists(dir_name)){
  dir.create(dir_name)
}

# loop through gene_ids to create a feature plot for each gene
for (gene_info in Zone_markers_t_gene_list ) {
  
  # create a feature plot for the gene
  feature_plot <- FeaturePlot(processed_LRP, features = gene_info[1], order = T, label = T, pt.size = 3) +
    scale_color_gradientn(colors = c('lightgray','yellow','red','darkred')) +
    ggtitle(paste0(gene_info[1], " - ", gene_info[2]))
  
  # Add the gene of interest to the metadata
  processed_LRP[["gene_of_interest"]] <- ifelse(processed_LRP[["RNA"]]@counts[gene_info[1],] > 0, gene_info[1], "Other")
  
  # create a dim plot for the gene
  dim_plot <- DimPlot(processed_LRP, group.by = "gene_of_interest") +
    ggtitle(paste0(gene_info[1], " - ", gene_info[2]))
  
  # Combine the plots into a single grob
  combined_plot <- gridExtra::arrangeGrob(feature_plot, dim_plot, ncol=2)
  
  # Specify the name of the PNG file
  png_file_name <- file.path(dir_name, paste0("Pericycle_phloem_pole_genes_", gene_info[1], "_", gene_info[2], ".png"))
  
  # Save the combined plot as a PNG file, specifying a larger width and height
  ggsave(filename = png_file_name, plot = combined_plot, width = 15, height = 5)
  
  # Remove the gene of interest from the metadata to prevent confusion in the next iteration
  processed_LRP[["gene_of_interest"]] <- NULL
}