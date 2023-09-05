library(abind)
library(conflicted)
library(cowplot)
library(dittoSeq)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(grid)
library(gridExtra)
library(magrittr)
library(patchwork)
library(png)
library(scales)
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(shiny)
library(stringr)
library(tidyr)
library(tidyverse)
library(shinyjs)
library(rsconnect)
library(cachem)

# Load the saved cc object
load("labeled_LRP.RData")

# Convert the labeled image to a data frame
cc_df <- as.data.frame(as.table(cc))

# Load processed LRP dataset
processed_LRP <- readRDS('LRP.processed2.rds')

# Assign cluster names to labels
cluster_names <- Idents(processed_LRP)
label_to_cluster <- data.frame(
  label = unique(cc_df$Freq), 
  cluster = rep(as.character(cluster_names), length.out = length(unique(cc_df$Freq)))
)

# Merge the label to cluster dataframe with the cc_df dataframe
cc_df <- merge(cc_df, label_to_cluster, by.x = "Freq", by.y = "label")

# Define the user interface
ui <- fluidPage(
  shinyjs::useShinyjs(),  
  titlePanel("Gene Expression Plot"),
  sidebarLayout(
    sidebarPanel(
      textInput("ATG_number", label = "Enter an ATG number:", value = "AT3G11260"),
      actionButton("plot_button", "Generate Plots"),
      textOutput("gene_error")  
    ),
    mainPanel(
      plotOutput("gene_expression_plot"),
      plotOutput("gene_plot")
    )
  )
)

# Define the server logic
server <- function(input, output) {
  output$gene_error <- renderText({ "" })  # Initialize the error message
  
  observeEvent(input$plot_button, {
    shinyjs::disable("plot_button")  # Disable the button
    
    withProgress(message = 'Calculating plots. Please wait...', {
      # Get the gene ID from the user input
      gene_id <- input$ATG_number
      
      # Check if gene is in cache
      cache <- gene_cache()
      if (gene_id %in% names(cache)) {
        output$gene_expression_plot <- renderPlot({ print(cache[[gene_id]]$gene_expression_plot) })
        output$gene_plot <- renderPlot({ print(cache[[gene_id]]$gene_plot) })
        shinyjs::enable("plot_button")  # Re-enable the button
        return() # Exit the function early
      }
      
      # Your existing code for checking if the gene exists in the dataset
      
      # Your existing code for calculating avg_gene_expression, merging with cc_df, etc.
      
      # Your existing code for creating p and feature_plot
      
      # Store the plots in the cache
      cache[[gene_id]] <- list(gene_expression_plot = p, gene_plot = feature_plot)
      gene_cache(cache)
      
      # Output the plots
      output$gene_expression_plot <- renderPlot({ print(p) })
      output$gene_plot <- renderPlot({ print(feature_plot) })
      
      shinyjs::enable("plot_button")  # Re-enable the button
    })
  })
}

# Run the application 
shinyApp(ui = ui, server = server)