library(EBImage)
library(tidyr)
library(dplyr)
library(ggplot2)
library(magick)


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