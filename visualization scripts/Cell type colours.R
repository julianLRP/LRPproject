library(EBImage)
library(tidyr)
library(dplyr)
library(ggplot2)
library(magick)
library(Seurat)
library(scales)
library(tidyverse)
library(reshape2)
library(png)
library(grid)

# Load the saved cc object
load("labeled_LRP.RData")

# Convert the labeled image to a data frame
cc_df <- as.data.frame(as.table(cc))


# Load the label-celltype mapping from CSV
label_celltype_df <- read.csv("Color_mapping_and_Identity.csv")

# Merge the label-celltype mapping with cc_df
cc_df <- merge(cc_df, label_celltype_df, by.x = "Freq", by.y = "Label")

# Reverse the levels of Var2
cc_df$Var2 <- factor(cc_df$Var2, levels = rev(levels(cc_df$Var2)))

# Define the color mapping for cell types
cell_type_colors <- c("0" = "white",  # background
                      "1" = "white",  # line segments
                      "Pericycle" = "#DAA060",  # beige
                      "Stele" = "#008000",      # green
                      "Endodermis" = "#FF0000", # red
                      "Epidermis" = "#FFA500", # orange
                      "QC" = "#FFFF00",         # yellow
                      "Cortex" = "#D8BFD8",     # light purple
                      "Protoxylem"="#5D3FD3",    #Iris
                      "LRC"="#7C4700",
                      "Columella"="#ADD8E5")      
                      

# Read the original image
img_for_plot <- readPNG('LRPstagesall3.png')

# Check if the image is grayscale
if (length(dim(img_for_plot)) == 2) {
  img_for_plot <- abind::abind(replicate(3, img_for_plot), array(1, dim = dim(img_for_plot)), along = 3)
} else if (dim(img_for_plot)[3] == 3) {
  img_for_plot <- abind::abind(img_for_plot, array(1, dim = dim(img_for_plot)[1:2]), along = 3)
} else if (dim(img_for_plot)[3] == 4) {
  img_for_plot[,,4] <- img_for_plot[,,4] * 1
}

# Create a raster graphic with transparency
g <- rasterGrob(img_for_plot, interpolate=TRUE)

# Set the X and Y location of the labels
labels_x <- c(185.0162, 446.6079, 753.2363, 1064.4975, 1384.5490, 1709.9236, 2013.7494, 2350.4541, 2777.3491)
labels_y <- 1400

# Create the ggplot for cell types
p_celltypes <- ggplot(cc_df, aes(x = Var1, y = Var2, fill = Cell_type
)) +  
  geom_raster(interpolate = TRUE) +
  scale_fill_manual(values = cell_type_colors, na.value = "white") +
  theme_void() +
  theme(plot.background = element_rect(fill = "white"),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        legend.background = element_rect(fill = "white", color = NA),
        legend.key = element_rect(fill = "white", color = NA),
        legend.box = element_blank(),
        legend.box.background = element_rect(colour = "white", linewidth = 0.2),
        legend.box.margin = margin(6, 6, 6, 6)) +
  coord_fixed() +
  annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
  labs(fill = "Cell Type") +
  annotate("text", 
           x = labels_x,
           y = rep(labels_y, length(labels_x)),
           label = c("pre-branch stage", "LRP 0", "LRP 1", "LRP 2", "LRP 3", "LRP 4", "LRP 5", "LRP 6", "LRP 7"), 
           size = 4, color = "black", hjust = 0.5, vjust = 1.5)

# Print the plot
print(p_celltypes)