**README for ColorSegment: A User-Customizable Image Labeling Application**
Welcome to ColorSegment, a robust image labeling tool designed to facilitate the process of color-based segmentation.

**Introduction**
ColorSegment enables you to assign custom color labels to the areas within your images. This intuitive application supports an unlimited number of labels, allowing you to handle complex images with diverse segmentation needs. Whether you're managing five labels or fifty, ColorSegment can accommodate you.

How To Use
The ease-of-use of ColorSegment will empower you to upload and label as many images as you'd like. you can manually color the segments by simply clicking on your picture.

Follow these steps to get started:

Upload your grayscale image in the application: ColorSegment supports only PNG

Interact with the user interface (UI): The UI is designed to be intuitive and user-friendly. Here's a quick look at what you can expect:

library(shiny)
library(png)

ui <- fluidPage(
  fileInput("upload", "Upload your image", accept = c('image/png', 'image/jpeg')),
  imageOutput("img"),
  colourInput("color", "Choose color"),
  textInput("label", "Enter label for this color"),
  actionButton("save", "Save label"),
  tableOutput("labels"),
  actionButton("save_image", "Save labeled image")
)

Upload, Choose Color, and Label: Use the fileInput to upload your image, the colourInput to choose your label color, and textInput to assign a name to the chosen color.

Save your label: Once done, use the save action button to save your label. You can view all labels using the labels table output.

Save your work: After all segments have been colored and labeled as desired, save your final, labeled image by clicking the save_image action button.
