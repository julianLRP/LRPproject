# image_processing.R

library(EBImage)

# Define a function to process an image
process_image <- function(img_path) {
  # Read the image
  img <- readImage(img_path)
  
  # Convert to grayscale if necessary
  if (length(dim(img)) > 2 && dim(img)[3] > 1) {
    img <- channel(img, 'gray')
  }
  
  # Normalize the image
  img <- normalize(img)
  
  # Compute Otsu's threshold on the grayscale image
  threshold <- otsu(img)
  
  # Apply the threshold to the image
  img <- img > threshold
  
  # Create a distance map of the thresholded image
  img <- distmap(img)
  
  # Perform a watershed transformation on the distance map
  img <- watershed(img)
  
  # Transform watershed image into binary image
  img <- ifelse(img > 1, 1, 0)
  
  # Label the watershed image
  img <- bwlabel(img)
  
  return(img)
}