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

# Define the user interface
ui <- fluidPage(
  shinyjs::useShinyjs(),
  titlePanel("Gene Expression Plot"),
  sidebarLayout(
    sidebarPanel(
      textInput("gene_input", label = "Enter a gene ID:", value = "AT5G26930"),
      actionButton("plot_button", "Generate Plots"),
      textOutput("gene_error")
    ),
    mainPanel(
      plotOutput("gene_expression_plot")
    )
  )
)

# Define the server logic
server <- function(input, output) {
  
  # Create a reactive value to store the gene data
  reactive_gene_data <- reactiveValues(data = NULL)
  
  observeEvent(input$plot_button, {
    shinyjs::disable("plot_button")  # Disable the button during calculation
    
    withProgress(message = 'Calculating plots. Please wait...', {
      # Load necessary data
      load("labeled_LRP.RData")
      processed_LRP <- readRDS('LRP.processed2.rds')
      mapping_df <- read.csv("Color_mapping_and_Identity.csv")
      
      # Extract earlymappingbycelltype values and merge with cc_df
      label_to_earlymapping <- mapping_df[, c("Label", "earlymappingbycelltype", "LRP_stage")]
      colnames(label_to_earlymapping)[which(names(label_to_earlymapping) == "Label")] <- "Freq"
      cc_df <- merge(as.data.frame(as.table(cc)), label_to_earlymapping, by = "Freq")
      
      # Calculate the average gene expression for each (earlymappingbycelltype)
      gene_of_interest <- input$gene_input  # Get the gene ID from the user input
      avg_gene_expression <- AverageExpression(object = processed_LRP, assays = "RNA", features = gene_of_interest, group.by = "earlymappingbycelltype")
      
      # Convert avg_gene_expression list into a data frame and reshape
      avg_gene_expression_df <- as.data.frame(avg_gene_expression) %>%
        gather(key = "earlymappingbycelltype", value = gene_of_interest) %>%
        mutate(earlymappingbycelltype = gsub("RNA.", "", earlymappingbycelltype))
      
      # Merge with cc_df based on "earlymappingbycelltype".
      cc_df <- merge(cc_df, avg_gene_expression_df, by = c("earlymappingbycelltype"), all.x = TRUE)
      
      # Set the gene_of_interest value of the pixels with label 0 or 1 to NA
      cc_df$gene_of_interest[cc_df$Freq %in% c(0, 1)] <- NA
      
      # Reverse the levels of Var2
      cc_df$Var2 <- factor(cc_df$Var2, levels = rev(levels(cc_df$Var2)))
      
      # Read the original image
      img_for_plot <- readPNG('LRPstagesall3.png')
      
      # Check if the image is grayscale or RGB, and adjust accordingly
      if (length(dim(img_for_plot)) == 2) { 
        # The image is grayscale. We need to convert it to an RGB image.
        img_for_plot <- abind::abind(replicate(3, img_for_plot), along = 3)
      } else if (dim(img_for_plot)[3] == 3) { 
        # The image is RGB. We need to add an alpha channel.
        img_for_plot <- abind::abind(img_for_plot, array(255, dim = c(dim(img_for_plot)[1:2], 1)), along = 3)
      }
      
      # If the image is already RGBA (with transparency channel), no changes are necessary.
      
      # Optionally, if you need to set a specific level of transparency, you can do so by modifying the alpha channel.
      # For example, to make the image 50% transparent, you could do:
      # img_for_plot[,,4] <- img_for_plot[,,4] * 0.5
      
      # Create a raster graphic object with the image, to be used in ggplot2
      g <- rasterGrob(img_for_plot, interpolate = TRUE)
      
      # Create a raster graphic with transparency
      g <- rasterGrob(img_for_plot, interpolate=TRUE)
      
      # Create the ggplot
      p <- ggplot(cc_df, aes(x = Var1, y = Var2, fill = gene_of_interest)) +
        geom_raster(interpolate = TRUE) +
        scale_fill_gradientn(colors = c('lightgray', 'yellow', 'red', 'darkred'), na.value = "white") +
        theme_void() +
        theme(
          plot.background = element_rect(fill = "white"),
          panel.border = element_blank(),  # Remove the plot box
          panel.grid = element_blank(),  # Remove grid lines
          legend.background = element_rect(fill = "white", color = NA),
          legend.key = element_rect(fill = "white", color = NA), 
          legend.box = element_blank(),  # Remove the legend box
          legend.box.background = element_rect(colour = "white", linewidth = 0.2),
          legend.box.margin = margin(6, 6, 6, 6)
        ) +
        coord_fixed() +
        annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +  # Add the raster image
        labs(fill = gene_of_interest) +  # Change the legend title
        annotate("text", 
                 x = labels_x,  # Place labels at the specific x coordinates
                 y = rep(labels_y, length(labels_x)),  # Place all labels at the same y-coordinate
                 label = c("LRP pre-branch stage", "LRP 0", "LRP 1", "LRP 2", "LRP 3", "LRP 4", "LRP 5", "LRP 6", "LRP 7"), 
                 size = 4, color = "black", hjust = 0.5, vjust = 1.5)
        
        # Store the final plot in the reactive value
        reactive_gene_data$data <- p
    })
    
    shinyjs::enable("plot_button")  # Re-enable the button after plot generation
  })
  
  # Render the plot
  output$gene_expression_plot <- renderPlot({
    if (!is.null(reactive_gene_data$data)) {
      print(reactive_gene_data$data)
    }
  })
}

# Run the application 
shinyApp(ui = ui, server = server)