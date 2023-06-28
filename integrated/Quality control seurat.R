#Load RDS
integerated_data <- readRDS("Integrated.processed.rds")

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
DimPlot(processed_LRP, label = TRUE)

if (!dir.exists("Elbowplots")) {
  dir.create("Elbowplots")
}

# Create the Elbow plot object
Elbow_plot <- ElbowPlot(integrated_data)

# Save the plot as a PNG image in the "Elbowplots" directory
ggsave(filename = "Elbowplots/elbow_plot.png", plot = Elbow_plot, dpi = 300)