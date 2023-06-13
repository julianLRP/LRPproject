# helpers.R

# Define a helper function to check if a file is an image
is_image_file <- function(filename) {
  # Get the file extension
  ext <- tools::file_ext(filename)
  
  # Check if the extension is a common image file extension
  return(tolower(ext) %in% c('jpg', 'jpeg', 'png', 'bmp', 'tif', 'tiff'))
}