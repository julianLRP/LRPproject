# Load required packages
library(shiny)
library(EBImage)

# Load custom functions and modules
source("image_processing.R")
source("user_interaction.R")

# Define the Shiny app
shinyApp(
  ui = ui, 
  server = server
)