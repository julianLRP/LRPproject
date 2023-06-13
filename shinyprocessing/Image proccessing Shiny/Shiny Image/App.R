#App.R
# Load required libraries
library(shiny)
library(shinyWidgets)

# Source the scripts
source("global.R")
source("helpers.R")
source("image_processing_shiny.R")
source("user_interaction.R")
source("server.R")

# Run the app
shinyApp(ui = ui, server = server)






