library(EBImage)
library(tidyr)
library(dplyr)
library(ggplot2)
library(magick)

# Read the image
img <- readImage("LRPstagesall.png")

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

# Define labels that should be considered as background
background_labels <- c(2,313, 0,322, 1)

# Set the selected labels as background
for(label in background_labels) {
  cc[cc == label] <- 0
}

writeImage(cc, "labeled_LRP.png")

save(cc, file = "labeled_LRP.RData")

display( colorLabels(cc))

#the following plot is to identify which label integer belongs to which segment via visualization.
# Convert the labeled image to a data frame
cc_df <- as.data.frame(as.table(cc))

# Compute the centroid of each label
centroids <- aggregate(cbind(Var1, Var2) ~ Freq, data = cc_df, FUN = function(x) mean(as.numeric(x)))

# Define your plot and assign it to 'p'
p <- ggplot(cc_df, aes(x = Var1, y = as.numeric(as.character(Var2)))) +  # removed fill = Freq
  geom_text(data = centroids, aes(x = Var1, y = as.numeric(as.character(Var2)), label = Freq), size = 0.5, color = "red") +
  theme_minimal() +
  theme(
    axis.text = element_blank(), 
    axis.title = element_blank(), 
    axis.ticks = element_blank(),
    panel.background = element_rect(fill = "white"),  # set background to white
    panel.grid = element_blank()  # remove grid lines
  ) +
  coord_fixed() +
  scale_y_reverse()  # Reverse the y-axis

# Save the plot
ggsave(filename = "labeled_image.png", plot = p, dpi = 500)
