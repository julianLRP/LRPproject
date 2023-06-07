# Load required packages
library(shiny)
library(EBImage)
library(png)
# Load custom functions and modules
source("image_processing.R")
source("user_interaction.R")

# Define the Shiny app
shinyApp(
  ui = ui <- fluidPage(
    fileInput("upload", "Upload your image", accept = c('image/png', 'image/jpeg')),
    imageOutput("img"),
    colourInput("color", "Choose color"),
    textInput("label", "Enter label for this color"),
    actionButton("save", "Save label"),
    tableOutput("labels"),
    actionButton("save_image", "Save labeled image")
    , 
  server = server
))
