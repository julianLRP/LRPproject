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

# Read the image
img <- readImage("LRPstagesall3.png")

# Check if the image is grayscale
if (length(dim(img)) == 2 || dim(img)[3] == 1) {
  # The image is grayscale. You can continue to thresholding.
  print("The image is already grayscale.")
} else {
  # The image is color. Convert it to grayscale.
  print("Converting the image to grayscale.")
  img <- channel(img, 'gray')
}

# Normalize the image
bw_img <- normalize(img)

# Compute Otsu's threshold on the grayscale image
threshold <- otsu(bw_img)

# Apply the threshold to the image
thresh_img <- bw_img > threshold

# Create a distance map of the thresholded image
dist_map <- distmap(thresh_img)

# Perform a watershed transformation on the distance map
ws <- watershed(dist_map)

# Transform watershed image into binary image
ws_binary <- ifelse(ws > 1, 1, 0)

# Label the watershed image
cc <- bwlabel(ws_binary)

# Save the labeled image
save(cc, file = "labeled_LRP.RData")

# Load the saved cc object
load("labeled_LRP.RData")

# Convert the labeled image to a data frame
cc_df <- as.data.frame(as.table(cc))

# Load processed LRP dataset
processed_LRP <- readRDS('Integrated.processed2.rds')

# Assign cluster names to labels randomly
cluster_names <- Idents(processed_LRP)
label_to_cluster <- data.frame(
  label = unique(cc_df$Freq), 
  cluster = rep(as.character(cluster_names), length.out = length(unique(cc_df$Freq)))
)

#Manual associations between labels and clusters
label_to_cluster <- data.frame(
  label = c("label1", "label2", "label3", ...),  # replace with your labels
  cluster = c("cluster1", "cluster2", "cluster3", ...)  # replace with associated clusters
)

# Merge the label to cluster dataframe with the cc_df dataframe
cc_df <- merge(cc_df, label_to_cluster, by.x = "Freq", by.y = "label")

# Calculate the average gene expression for each pixel
gene_of_interest <- "AT3G11260"
avg_gene_expression <- AverageExpression(object = processed_LRP, assays = "RNA", features = gene_of_interest)

# convert avg_gene_expression list into a data frame
avg_gene_expression <- as.data.frame(avg_gene_expression)

#reshape this data frame from wide format to long format, making sure to have one row per cluster
avg_gene_expression <- avg_gene_expression %>% 
  gather(key = "cluster", value = gene_of_interest)

#remove the "RNA." prefix from the cluster column
avg_gene_expression$cluster <- gsub("RNA.", "", avg_gene_expression$cluster)

# Merge the avg_gene_expression with cc_df
cc_df <- merge(cc_df, avg_gene_expression, by = "cluster")

# Scale the gene_expression column to between 0 and 1
cc_df$gene_expression <- scales::rescale(cc_df$gene_expression, to = c(0, 1))

# Set the gene_of_interest value of the pixels with label 0 or 1 to NA this will set all values of the background and the lines to NA.
cc_df$gene_of_interest[cc_df$Freq %in% c(0, 1)] <- NA

# Scale the gene_expression column to between 0 and 1
cc_df$gene_expression <- scales::rescale(cc_df$gene_expression, to = c(0, 1))

# Reverse the levels of Var2
cc_df$Var2 <- factor(cc_df$Var2, levels = rev(levels(cc_df$Var2)))

# Read the original image
img_for_plot <- readPNG('LRPstagesall3.png')

##will execute the overlaying of the original image regardless of gray scale or RGBA image
# Check if the image is grayscale
if (length(dim(img_for_plot)) == 2) {
  # Convert grayscale image to RGBA image
  img_for_plot <- abind::abind(replicate(3, img_for_plot), array(1, dim = dim(img_for_plot)), along = 3)
} else if (dim(img_for_plot)[3] == 3) {
  # Add an alpha channel to RGB image
  img_for_plot <- abind::abind(img_for_plot, array(1, dim = dim(img_for_plot)[1:2]), along = 3)
} else if (dim(img_for_plot)[3] == 4) {
  # Modify the existing alpha channel of RGBA image
  img_for_plot[,,4] <- img_for_plot[,,4] * 1  # Adjust transparency level here
}

# Create a raster graphic with transparency
g <- rasterGrob(img_for_plot, interpolate=TRUE)

#Set the X an Y location of the labels
labels_x <- c(132.0162, 446.6079, 753.2363, 1064.4975, 1384.5490, 1709.9236, 2013.7494, 2350.4541, 2777.3491)
labels_y <- 1400#keep between Y axis of original image to prevent annotation custom from being misaligned with plotted label pixels

# Create the ggplot
p <- ggplot(cc_df, aes(x = Var1, y = Var2, fill = gene_of_interest)) +
  geom_raster(interpolate = TRUE) +
  scale_fill_gradientn(colors = c('lightgray', 'yellow', 'red', 'darkred'), na.value = "white") +
  theme_void() +
  theme(plot.background = element_rect(fill = "white"),
        panel.border = element_blank(),  # Remove the plot box here
        panel.grid = element_blank(),  # Remove grid lines
        legend.background = element_rect(fill = "white", color = NA),
        legend.key = element_rect(fill = "white", color = NA), 
        legend.box = element_blank(),  # Remove the legend box here
        legend.box.background = element_rect(colour = "white", linewidth = 0.2),
        legend.box.margin = margin(6, 6, 6, 6)) +
  coord_fixed() +
  annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +  # Add this line
  labs(fill = gene_of_interest) +  # Change the legend title here
  annotate("text", 
           x = labels_x,  # Place labels at the specific x coordinates
           y = rep(labels_y, length(labels_x)),  # Place all labels at the same y-coordinate
           label = c("LRP pre-branch stage", "LRP 0", "LRP 1", "LRP 2", "LRP 3", "LRP 4", "LRP 5", "LRP 6", "LRP 7"), 
           size = 4, color = "black", hjust = 0.5, vjust = 1.5)

# Print the plot
print(p)
