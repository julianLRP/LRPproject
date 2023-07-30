# Load the required libraries
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

# Define the user interface
ui <- fluidPage(
  titlePanel("Gene Expression Plot"),
  sidebarLayout(
    sidebarPanel(
      textInput("ATG_number", label = "Enter an ATG number:", value = "AT3G11260"),
      actionButton("plot_button", "Generate Plots")
    ),
    mainPanel(
      plotOutput("gene_expression_plot"),
      plotOutput("gene_plot")
    )
  )
)

# Define the server logic
server <- function(input, output) {
  observeEvent(input$plot_button, {
    # Get the gene ID from the user input
    gene_id <- input$ATG_number
    
    # Load the saved cc object
    load("labeled_LRP.RData")
    
    # Convert the labeled image to a data frame
    cc_df <- as.data.frame(as.table(cc))
    
    # Load processed LRP dataset
    processed_LRP <- readRDS('Integrated.processed2.rds')
    
    # Assign cluster names to labels
    cluster_names <- Idents(processed_LRP)
    label_to_cluster <- data.frame(
      label = unique(cc_df$Freq), 
      cluster = rep(as.character(cluster_names), length.out = length(unique(cc_df$Freq)))
    )
    
    # Merge the label to cluster dataframe with the cc_df dataframe
    cc_df <- merge(cc_df, label_to_cluster, by.x = "Freq", by.y = "label")
    
    # Calculate the average gene expression for each pixel
    gene_of_interest <- gene_id
    avg_gene_expression <- AverageExpression(object = processed_LRP, assays = "RNA", features = gene_of_interest)
    
    # convert avg_gene_expression list into a data frame
    avg_gene_expression <- as.data.frame(avg_gene_expression)
    
    #reshape this data frame from wide format to long format, making sure we have one row per cluster
    avg_gene_expression <- avg_gene_expression %>% 
      gather(key = "cluster", value = gene_of_interest)
    
    #remove the "RNA." prefix from the cluster column
    avg_gene_expression$cluster <- gsub("RNA.", "", avg_gene_expression$cluster)
    
    # Merge the avg_gene_expression with cc_df
    cc_df <- merge(cc_df, avg_gene_expression, by = "cluster")
    
    # Scale the gene_expression column to between 0 and 1
    cc_df$gene_expression <- scales::rescale(cc_df$gene_expression, to = c(0, 1))
    
    # Set the gene_of_interest value of the pixels with label 0 or 1 to NA
    cc_df$gene_of_interest[cc_df$Freq %in% c(0, 1)] <- NA
    
    # Scale the gene_expression column to between 0 and 1
    cc_df$gene_expression <- scales::rescale(cc_df$gene_expression, to = c(0, 1))
    
    # Reverse the levels of Var2
    cc_df$Var2 <- factor(cc_df$Var2, levels = rev(levels(cc_df$Var2)))
    
    # Read the original image
    img_for_plot <- readPNG('LRPstagesall3.png')
    
    # Create a raster graphic with transparency
    g <- rasterGrob(img_for_plot, interpolate=TRUE)
    
    #Set the X an Y location of the labels
    labels_x <- c(132.0162, 446.6079, 753.2363, 1064.4975, 1384.5490, 1709.9236, 2013.7494, 2350.4541, 2777.3491)
    labels_y <- 1400
    
    # Create the ggplot
    p <- ggplot(cc_df, aes(x = Var1, y = Var2, fill = gene_of_interest)) +
      geom_raster(interpolate = TRUE) +
      scale_fill_gradientn(colors = c('lightgray', 'yellow', 'red', 'darkred'), na.value = "white") +
      theme_void() +
      theme(plot.background = element_rect(fill = "white"),
            panel.border = element_blank(),  # Remove the plot box here
            panel.grid = element_blank(),  # Remove grid lines
            legend.background = element_rect(fill = "white", color = NA),
            legend.key = element_rect(fill = "white", color = NA), 
            legend.box = element_blank(),  # Remove the legend box here
            legend.box.background = element_rect(colour = "white", linewidth = 0.2),
            legend.box.margin = margin(6, 6, 6, 6)) +
      coord_fixed() +
      annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +  # Add this line
      labs(fill = gene_of_interest) +  # Change the legend title here
      annotate("text", 
               x = labels_x,  # Place labels at the specific x coordinates
               y = rep(labels_y, length(labels_x)),  # Place all labels at the same y-coordinate
               label = c("LRP pre-branch stage", "LRP 0", "LRP 1", "LRP 2", "LRP 3", "LRP 4", "LRP 5", "LRP 6", "LRP 7"), 
               size = 4, color = "black", hjust = 0.5, vjust = 1.5)
    
    # Output the plot
    output$gene_expression_plot <- renderPlot({print(p)})
    
    # Read the CSV file
    gene_info_data <- read.csv("gene_info.csv", stringsAsFactors = FALSE)
    
    # Convert the data frame to a named list
    gene_descriptions <- setNames(gene_info_data$gene_description, gene_info_data$gene_id)
    
    # Define the function to get the gene description
    get_gene_description <- function(gene_id) {
      return(gene_descriptions[[as.character(gene_id)]])
    }
    
    # Check if the gene exists in the dataset
    if(!gene_id %in% rownames(processed_LRP[["RNA"]]@counts)) {
      print(paste("Skipping gene_id", gene_id, "as it's not present in the dataset."))  # Print message to the console
      return()  # Exit the function early
    }
    
    # Check if the gene exists in gene_info.csv
    if(!gene_id %in% gene_info_data$gene_id) {
      print(paste("Skipping gene_id", gene_id, "as it's not present in the gene_info.csv."))  # Print message to the console
      return()  # Exit the function early
    }
    
    # Get the gene description
    gene_description <- get_gene_description(gene_id)
    
    # Get the gene description
    gene_description <- get_gene_description(gene_id)
    
    # Check if the gene exists in the dataset
    if(!gene_id %in% rownames(processed_LRP[["RNA"]]@counts)) {
      print(paste("Skipping gene_id", gene_id, "as it's not present in the dataset."))  # Print message to the console
    }
    
    # Get UMAP data and convert it to a data frame
    umap_data <- as.data.frame(Embeddings(processed_LRP, "umap"))
    names(umap_data) <- c("UMAP1", "UMAP2")  # Explicitly name the columns
    
    # Add cell names and cluster IDs to the data frame
    umap_data$cell <- rownames(umap_data)
    umap_data$cluster <- Idents(processed_LRP)
    
    # Calculate the centroids of the clusters
    cluster_centroids <- umap_data %>%
      group_by(cluster) %>%
      summarize(UMAP1 = mean(UMAP1), UMAP2 = mean(UMAP2))
    
    # Fetch the gene expression data
    gene_data <- as.data.frame(GetAssayData(processed_LRP, slot = "data")[gene_id, ])
    names(gene_data) <- c("Expression")  # Assuming that gene_data has only one column
    gene_data$cell <- rownames(gene_data)
    
    # Combine UMAP and gene expression data
    combined_data <- full_join(umap_data, gene_data, by = "cell")
    
    # Create a feature plot for the gene
    feature_plot <- ggplot(combined_data, aes(x = UMAP1, y = UMAP2)) +
      geom_point(aes(color = Expression), size = 1.5) +
      scale_color_gradientn(colors = c('lightgray','yellow','red','darkred')) +
      geom_text_repel(data = cluster_centroids, aes(label = cluster), size = 3, box.padding = unit(0.35, "lines")) +
      ggtitle(paste0(gene_id)) +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            panel.background = element_blank())
    
    # Add the gene of interest to the metadata
    processed_LRP[["gene_of_interest"]] <- ifelse(processed_LRP[["RNA"]]@counts[gene_id,] > 0, gene_id, "Other")
    
    # Create a dim plot for the gene
    color_vector <- setNames(c(rgb(1,0,0,alpha=0.5), rgb(0.5,0.5,0.5,alpha=0.5)), c(gene_id, "Other"))
    dim_plot <- DimPlot(processed_LRP, group.by = "gene_of_interest", pt.size = 3, cols = color_vector) +
      ggtitle(paste0(gene_id))
    
    # Create a description text box
    wrapped_gene_description <- str_wrap(gene_description, width = 80)  # Adjust the width as needed
    description_lines <- str_split(wrapped_gene_description, "\n", simplify = FALSE)[[1]]
    line_count <- length(description_lines)
    
    # Create a data frame with the lines and their corresponding y-positions
    y_increment <- 1 / (line_count + 1)  # Adjust this value to increase or decrease the space between lines
    description_df <- data.frame(
      x = rep(0, line_count),
      y = seq(from = 1 - y_increment, to = 0, length.out = line_count),
      text = description_lines
    )
    
    # Calculate font size based on the number of lines
    if (line_count <= 10) {
      font_size <- 5
    } else if (line_count <= 20) {
      font_size <- 4
    } else {
      font_size <- 3
    }
    
    description_text <- ggplot(data = description_df) +
      geom_text(aes(x = x, y = y, label = text), hjust = 0, vjust = 0, size = font_size, family = "mono") +  # Use the calculated font_size
      theme_void() +
      theme(plot.margin = margin(t = 5, r = 5, b = 5, l = 5, unit = "pt"),
            text = element_text(family = "mono")) +
      coord_cartesian(clip = 'off', xlim = c(0, 1), ylim = c(0, 1))
    
    # Combine the plots and the description text box
    combined_plot <- (feature_plot / dim_plot / description_text) +
      patchwork::plot_layout(ncol = 1, heights = c(2, 2, 1))
    
    # Remove the gene of interest from the metadata to prevent confusion in the next iteration
    processed_LRP[["gene_of_interest"]] <- NULL
    
    # Output the plot
    output$gene_plot <- renderPlot({print(combined_plot)})
  })
}

# Run the application 
shinyApp(ui = ui, server = server)