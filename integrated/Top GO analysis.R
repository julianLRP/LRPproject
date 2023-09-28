##TopGO in combination with org.At.tair.db to find DEG's.

# Load required libraries
library(topGO)
library(org.At.tair.db)

# Load your Seurat object and get all genes
processed_LRP <- readRDS("Integrated.processed.rds")
gene_identifiers <- rownames(processed_LRP[["RNA"]]@counts)

# Prepare a dataframe for the gene mapping
gene_mapping <- data.frame(identifier = gene_identifiers,
                           gene_name = mapIds(org.At.tair.db, 
                                              keys = gene_identifiers, 
                                              keytype = "TAIR", 
                                              column = "SYMBOL"))

# Get all genes as a numeric vector, uninteresting by default
all_genes <- rep(0, length(gene_mapping$identifier))
names(all_genes) <- gene_mapping$identifier

# Define the cluster names (based on your CSV files)
clusters <- 0:36  # This will create a vector with numbers from 0 to 36

# Loop through each cluster
for (i in clusters) {
  
  # Read the marker genes CSV file
  markers_cluster <- read.csv(paste0(i, "_markers.csv"), stringsAsFactors = FALSE)
  
  # Define the interesting genes
  interesting_genes <- markers_cluster$identifier
  
  # Set the interesting genes to 1
  all_genes[interesting_genes] <- 1
  
  # Create a function for selecting the interesting genes
  geneSelFunc <- function(all_genes) { all_genes %in% interesting_genes }
  
  # Create a topGOdata object
  GOdata <- new("topGOdata",
                ontology = "BP",                  # Change to "CC" or "MF" for Cellular Component or Molecular Function
                allGenes = all_genes,
                geneSel = geneSelFunc,
                nodeSize = 10,
                annot = annFUN.org,
                mapping="org.At.tair.db",
                ID = "TAIR")
  
  # Run the GO enrichment analysis
  resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  
  # Get the top results
  topTableFisher <- GenTable(GOdata, classic = resultFisher, orderBy = "fisher", ranksOf = "classic", topNodes = 10)
  
  # Write the top table to a CSV file
  write.csv(topTableFisher, file = paste0("cluster", i, "_topGO.csv"), row.names = FALSE)
  
  # Reset all_genes to 0
  all_genes[interesting_genes] <- 0
}
