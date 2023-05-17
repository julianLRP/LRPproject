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
  library(SingleCellExperiment)


  # Load the saved sce object
  sce <- readRDS("LRP_with_pseudotime.rds")
  
  # Get the gene expression data
  counts <- assays(sce)$counts
  
  # Get the pseudotime values for each cell
  pseudotime <- colData(sce)$pseudotime
  
  # Get the gene expression values for the gene with ID "AT1G31950"
  gene_counts <- counts["AT1G31950",]
  
  # Create a data frame for plotting
  df <- data.frame(Pseudotime = pseudotime, Expression = gene_counts)
  
  # Plot the gene expression over pseudotime
  ggplot(df, aes(x = Pseudotime, y = Expression)) +
    geom_point() +
    theme_bw() +
    xlab("Pseudotime") +
    ylab("Expression of AT1G31950") +
    ggtitle("Expression of AT1G31950 over Pseudotime")